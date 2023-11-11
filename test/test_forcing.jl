using DEBmicroTrait, Test
using JLD2
using Plots, Dates
using OrdinaryDiffEq

# I/O
dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd();
avena_ex                = load(joinpath(dir, "files/avena/avena_exudation_2006_2022.jld2"));
avena_root              = load(joinpath(dir, "files/avena/avena_root_2006_2022.jld2"));
γ_root                  = 5e-5

# Plot
month_hours = 365.25/12*24
months = 12
nyears  = 2022-2006
tstop = month_hours*months*nyears
MONTHS = [Dates.LOCALES["english"].months... Dates.LOCALES["english"].months...]
xticks_month=(repeat((0:month_hours:tstop),nyears), repeat(MONTHS[1:months], nyears))
xticks_year = ((0:month_hours*12:month_hours*12*nyears), string.(collect(range(start=2006,step=1,stop=2022))))
plot(avena_ex["exudation_rate"][1:8760], xticks=xticks_month, xrotation=90)
plot(avena_ex["exudation_rate"][1:8760])

plot(avena_root["root_C"], xticks=xticks_year,  xrotation = 90)


# ODE
p_force = Forcing(avena_root["root_C"], γ_root, avena_ex["exudation_rate"])

function test!(du,u,p_force,t)
    spl_ex   = DEBmicroTrait.root_exudation(p_force)
    spl_root = DEBmicroTrait.root_decay(p_force)
    #
    du[1] = spl_ex(t)
    du[2] = spl_root(t)
end

tspan = [0,tstop]
u0 = [0.0, 0.0]
prob = ODEProblem(test!,u0,tspan,p)
sol = solve(prob, Tsit5());
plot(sol, idxs=2, xticks=xticks_year)





avena_TOC               = CSV.read(joinpath(dir, "files/avena/TOC_data_Hopland_Avena.csv"), DataFrame, missingstring="N/A")
avena_TOC_mM            = avena_TOC.TOC./12.011