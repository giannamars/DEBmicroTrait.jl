using CSV
using DataFrames
using JLD2               # JLD2 is the modern format; JLD works too
using DifferentialEquations
using Statistics
using ProgressMeter      # optional, nice progress bar
using Base.Threads       # for multi‑threading
using LinearAlgebra      # for @views, @inbounds
using DEBmicroTrait      # your custom DEB‑micro‑trait package
using Plots
gr()

"""
    run_isolate(id_isolate, id_monomer, id_polymer,
                assimilation, enzymes, maintenance,
                protein_synthesis, turnover;
                tspan=(0.0, 1000.0), n_save=1000)

Runs the DEB batch model for a single genome (`id_isolate`) and returns a
named tuple with the three medians (`r_median`, `Bio_median`, `BGE_median`).

All arguments after `turnover` are the dictionaries loaded from the JLD2
files.  The keyword arguments control the integration interval and the
number of saved points (default 1000 → matches your original script).
"""
function run_isolate(id_isolate::Int,
                     id_monomer::Int,
                     id_polymer::Int,
                     assimilation,
                     enzymes,
                     maintenance,
                     protein_synthesis,
                     turnover;
                     tspan = (0.0, 1000.0),
                     n_save = 1000)

    # ----------------------------------------------------------------
    # 1️⃣  Initialise the model for this isolate
    # ----------------------------------------------------------------
    p = DEBmicroTrait.init_batch_model(id_isolate,
                                       id_monomer,
                                       id_polymer,
                                       assimilation,
                                       enzymes,
                                       maintenance,
                                       protein_synthesis,
                                       turnover)

    # ----------------------------------------------------------------
    # 2️⃣  Build the initial condition vector `u0`
    # ----------------------------------------------------------------
    dim = p.setup_pars.dim
    n_polymers        = p.setup_pars.n_polymers
    n_monomers        = p.setup_pars.n_monomers
    n_microbes        = p.setup_pars.n_microbes
    n_enzymes         = p.setup_pars.n_enzymes
    u0 = zeros(dim)

    # The indices are a bit cryptic – they come from the original script.
    # I kept the same logic but used `@views` to avoid temporary copies.
    @views begin
        # Biomass pools (microbes) – 90 % in the “active” pool, 10 % in the “reserve”
        u0[1 + n_polymers + n_monomers : n_polymers + n_monomers + n_microbes] .=
            0.9 * initb["Bio0"][id_isolate]

        u0[1 + n_polymers + n_monomers + n_microbes : n_polymers + n_monomers + 2 * n_microbes] .=
            0.1 * initb["Bio0"][id_isolate]

        # Substrate concentrations (monomers) – constant 1.25 for all
        u0[1 + n_polymers : n_polymers + n_monomers] .= 1.25

        # Polymer concentrations – constant 10.0 for all
        u0[1 : n_polymers] .= 10.0

        # Enzyme pools – start at zero (already the default)
        # (kept for completeness)
        # u0[ … ] .= 0.0
    end

    # ----------------------------------------------------------------
    # 3️⃣  Define the ODE problem
    # ----------------------------------------------------------------
    prob = ODEProblem(DEBmicroTrait.batch_model!, u0, tspan, p)

    # Use a stiff solver (Rodas5 is a good default) and request exactly
    # `n_save` equally spaced output points.
    save_times = range(tspan[1], tspan[2]; length=n_save)

    sol = solve(prob,
                alg_hints=[:stiff];                     # stiff solver
                saveat = save_times,
                callback = termination_callback(),
                reltol = 1e-6,
                abstol = 1e-8)

    # ----------------------------------------------------------------
    # 4️⃣  Extract the three time‑series
    # ----------------------------------------------------------------
    ts = extract_time_series(sol, p)

    # ----------------------------------------------------------------
    # 5️⃣  Compute the medians (skip zeros / missing values)
    # ----------------------------------------------------------------
    r_med   = safe_median(ts.r)
    Bio_med = safe_median(ts.Bio)
    BGE_med = safe_median(ts.BGE)

    return (r_median = r_med,
            Bio_median = Bio_med,
            BGE_median = BGE_med)
end

# --------------------------------------------------------------------
# Callback that stops the integration when the first state variable
# (presumably the substrate) falls below 1e‑5.
# --------------------------------------------------------------------
function termination_callback()
    condition(u, t, integrator) = u[1] - 1e-5
    affect!(integrator) = terminate!(integrator)
    ContinuousCallback(condition, affect!)
end

# --------------------------------------------------------------------
# Compute the three quantities we need from a solved ODE.
# All arguments are vectors of the same length (the number of saved
# time points).  The function returns a NamedTuple with the three
# time‑series.
# --------------------------------------------------------------------
function extract_time_series(sol, p)
    # `sol.u` is a Vector{Vector{Float64}} – each element is the state
    # at a saved time point.
    n = length(sol.u)

    # Pre‑allocate output vectors
    r  = Vector{Float64}(undef, n)
    Bio = Vector{Float64}(undef, n)
    BR = Vector{Float64}(undef, n)
    BP = Vector{Float64}(undef, n)

    # Temporary derivative vector – reused for every call to batch_model!
    du = similar(sol.u[1])

    @inbounds for i in 1:n
        ui = sol.u[i]                     # state at time i
        # ---- growth rate -------------------------------------------------
        # `growth!` expects a 1‑element vector for the first argument,
        # the metabolism parameters, and the two substrate concentrations.
        # It returns a 1‑element vector; we take the first entry.
        r[i] = DEBmicroTrait.growth!(zeros(1), p.metabolism_pars,
                                     [ui[3]], [ui[4]])[1]

        # ---- biomass ----------------------------------------------------
        Bio[i] = ui[3] + ui[4]            # sum of the two biomass pools

        # ---- respiration & production ------------------------------------
        # `batch_model!` returns a vector of rates; the last entry is
        # respiration (BR) and entries 3 & 4 are the two production
        # fluxes that together give BP.
        rates = DEBmicroTrait.batch_model!(du, ui, p, 0.0)
        BR[i] = rates[end]               # respiration
        BP[i] = rates[3] + rates[4]       # production (both pools)
    end

    BGE = @. BP / (BP + BR)               # bacterial growth efficiency

    return (r=r, Bio=Bio, BGE=BGE)
end

# --------------------------------------------------------------------
# Compute the median of a vector while safely handling empty or all‑zero
# vectors.  Returns `missing` when no valid value exists.
# --------------------------------------------------------------------
function safe_median(v::AbstractVector{T}) where {T<:Real}
    # Remove zeros (or any other sentinel you prefer) and missing values
    filtered = filter(x -> x != 0 && !ismissing(x), v)
    return isempty(filtered) ? missing : median(filtered)
end

# --------------------------------------------------------------------
# 0️⃣  Paths & data loading
# --------------------------------------------------------------------
dir = get(ENV, "DEBSCRIPTS", pwd())

df_mags = CSV.read(joinpath(dir, "files/input/greenlon-H-mags2traits.csv"),
                  DataFrame; missingstring = "")

df_metabolites = CSV.read(joinpath(dir, "files/input/greenlon-H-monomers.csv"),
                         DataFrame; missingstring = "")

# Load the JLD2 dictionaries (they are simple `Dict{String,Any}` objects)
assimilation      = load(joinpath(dir, "files/output/mags_assimilation.jld"))
enzymes           = load(joinpath(dir, "files/output/mags_enzymes.jld"))
maintenance       = load(joinpath(dir, "files/output/mags_maintenance.jld"))
protein_synthesis = load(joinpath(dir, "files/output/mags_protein_synthesis.jld"))
turnover          = load(joinpath(dir, "files/output/mags_turnover.jld"))
initb             = load(joinpath(dir, "files/output/mags_batch_init.jld"))

# --------------------------------------------------------------------
# 1️⃣  Settings for the simulation
# --------------------------------------------------------------------
id_polymer = 2   # e.g. chitin
id_monomer = 2

n_isolates = nrow(df_mags)

# --------------------------------------------------------------------
# 2️⃣  Run all isolates – choose serial or threaded execution
# --------------------------------------------------------------------
results = Vector{NamedTuple}(undef, n_isolates)   # pre‑allocate

# Progress bar (nice visual feedback)
pbar = Progress(n_isolates; desc = "Simulating isolates", dt = 0.5)

# ---- Serial version (simpler, works on any Julia installation) ----
for i in 1:n_isolates
    results[i] = run_isolate(i, id_monomer, id_polymer,
                             assimilation, enzymes, maintenance,
                             protein_synthesis, turnover)
    next!(pbar)
end

# ---- Threaded version (uncomment to use all available threads) ----
# @threads for i in 1:n_isolates
#     results[i] = run_isolate(i, id_monomer, id_polymer,
#                              assimilation, enzymes, maintenance,
#                              protein_synthesis, turnover)
#     @sync @async next!(pbar)   # progress bar is thread‑safe via @async
# end

# --------------------------------------------------------------------
# 3️⃣  Assemble the results back into the original DataFrame
# --------------------------------------------------------------------
# Convert the vector of NamedTuples into a DataFrame
df_res = DataFrame(results)

# The original script normalised the biomass median by its own maximum.
# Guard against a zero‑maximum (which would give Inf/NaN).
max_bio = maximum(skipmissing(df_res.Bio_median))
if max_bio == 0 || isnan(max_bio)
    relabund = fill(missing, n_isolates)
else
    relabund = df_res.Bio_median ./ max_bio
end

# Add the new columns to `df_mags`
df_mags.r_median            = df_res.r_median
df_mags.relabund_median     = relabund
df_mags.bge_median          = df_res.BGE_median

# --------------------------------------------------------------------
# 4️⃣  Write the final table
# --------------------------------------------------------------------
out_path = joinpath(dir, "files/output/greenlon-H-chitin.csv")
CSV.write(out_path, df_mags)

println("\n✅  Finished! Results written to: $out_path")


using ForwardDiff
du = similar(u0)
J = ForwardDiff.jacobian(u -> begin
        DEBmicroTrait.batch_model!(u0, u, 0.0, p)   # note: `du` is mutated, we ignore it
        du
    end, u0)