using DifferentialEquations
using SciMLSensitivity
using Optimization, OptimizationOptimisers
using Flux
using Zygote
using Statistics
using DelimitedFiles
using Random
using Plots

# ==============================================================================
# 1. SETUP DATA
# ==============================================================================

# A. Load Z-Scores (1 x N)
# Assuming you have saved this from the previous step
z_scores = readdlm("genome_complexity_score.csv", ',', Float32)
z_scores = permutedims(z_scores) # Ensure shape (1, N)
n_genomes = size(z_scores, 2)

# B. Load Real qSIP Growth Rates (Vector of length N)
# For this demo, I generate synthetic data. REPLACE WITH YOUR CSV.
# qSIP rates usually range from 0.01 to 0.5 per day.
real_qsip_rates = Float32.(0.05 .+ 0.2 .* rand(n_genomes)) 

# C. Experimental Constants (CRITICAL)
# You must set these to match your actual experiment.
incubation_time = 24.0f0      # e.g., 24 hours
initial_substrate = 10.0f0    # e.g., 10 mg/L glucose added
initial_biomass = 0.1f0       # Assumed starting biomass (normalized)

t_span = (0.0f0, incubation_time)

println("Setup: $n_genomes genomes. Incubation: $incubation_time hours.")

# ==============================================================================
# 2. DEFINE THE MODEL
# ==============================================================================

# A. The ODE (Monod Kinetics)
function microbial_growth!(du, u, p, t)
    X, S = u
    # Parameters predicted by NN
    Vmax, K, m, Y = p
    
    # Physics
    S_safe = max(S, 0.0f0)
    uptake = Vmax * (S_safe / (S_safe + K))
    growth = (uptake * Y) - m
    
    du[1] = growth * X
    du[2] = -uptake * X
end

# B. The Neural Network
# Input: 1 (Z-Score) -> Output: 4 (DEB Params)
# We use a constrained architecture to ensure positive parameters
ann = Chain(
    Dense(1 => 16, tanh),
    Dense(16 => 16, tanh),
    # Softplus ensures Vmax, K, m, Y are strictly positive
    Dense(16 => 4, softplus) 
)

# Structure for Optimization
p_init, re = Flux.destructure(ann)

# ==============================================================================
# 3. PREDICTION & LOSS
# ==============================================================================

function predict_growth_rates(p_nn, z_batch)
    # 1. Decode Neural Net
    model = re(p_nn)
    bio_params = model(z_batch) # Shape (4, N)
    
    # 2. Define a helper function that solves for ONE genome index 'i'
    # This keeps the logic self-contained and avoids external mutation.
    function solve_single_genome(i)
        # Extract params for this specific genome
        p_ode = bio_params[:, i] # [Vmax, K, m, Y]
        
        # Initial Conditions (Create new vector every time, do not mutate old one)
        u0 = [initial_biomass, initial_substrate]
        
        # Define Problem
        prob = ODEProblem(microbial_growth!, u0, t_span, p_ode)
        
        # Solve
        sol = solve(prob, Tsit5(), save_everystep=false, sensealg=InterpolatingAdjoint())
        
        # Check success (Using simple if/else is differentiable)
        if sol.retcode == ReturnCode.Success
            X_final = sol[1, end]
            # Calculate Rate
            return (log(max(X_final, 1f-6)) - log(initial_biomass)) / incubation_time
        else
            return 0.0f0 # Fail penalty
        end
    end

    # 3. USE MAP INSTEAD OF LOOP/PUSH
    # Zygote can differentiate through 'map' perfectly.
    # We map the solver over every index 1..N
    predicted_rates = map(solve_single_genome, 1:size(z_batch, 2))
    
    # Ensure result is a concrete Float32 array
    return convert(Array{Float32}, predicted_rates)
end

function loss_function(p_nn, nothing)
    # Predict
    pred_rates = predict_growth_rates(p_nn, z_scores)
    
    # Compare with Observed qSIP rates (MSE)
    loss = mean(abs2, pred_rates .- real_qsip_rates)
    
    # Optional: Regularization to prevent extreme parameter values
    # e.g. penalize very large Ks or Vmaxs
    # reg = sum(abs2, p_nn) * 0.001 
    
    return loss
end

# ==============================================================================
# 4. TRAIN
# ==============================================================================

optf = OptimizationFunction(loss_function, Optimization.AutoZygote())
optprob = OptimizationProblem(optf, p_init)

println("Training Neural ODE... (Predicting Parameters from Complexity)")

callback = function (p, l)
    println("Loss: $l")
    return false
end

# Use Adam with a small learning rate (ODE gradients can be unstable)
result = Optimization.solve(optprob, OptimizationOptimisers.Adam(0.005), callback=callback, maxiters=100)

# ==============================================================================
# 5. RESULT EXTRACTION
# ==============================================================================

final_model = re(result.u)

# Sort Z-scores for plotting scaling laws
sorted_indices = sortperm(vec(z_scores))
z_sorted = z_scores[:, sorted_indices]
params_sorted = final_model(z_sorted) # Shape (4, N)

# Plot the learned biological laws
p1 = plot(vec(z_sorted), vec(params_sorted[1, :]), title="Vmax (Rate)", ylabel="/h", label=:none)
p2 = plot(vec(z_sorted), vec(params_sorted[2, :]), title="K (Half Sat)", ylabel="mg/L", label=:none)
p3 = plot(vec(z_sorted), vec(params_sorted[3, :]), title="m (Maintenance)", ylabel="/h", label=:none)
p4 = plot(vec(z_sorted), vec(params_sorted[4, :]), title="Y (Yield)", ylabel="g/g", label=:none)

display(plot(p1, p2, p3, p4, layout=(2,2), size=(800, 600)))

println("Done! The plots show how DEB parameters scale with Complexity.")


using Plots

# 1. Recover the Trained Model
# result.u contains the final optimized weights of the Neural Net
final_weights = result.u
trained_model = re(final_weights)

# 2. Sort Z-scores for clean plotting
# We want to plot from Smallest Genome (Low Z) to Largest (High Z)
perm = sortperm(vec(z_scores))
z_sorted = z_scores[:, perm] # Shape (1, N)

# 3. Predict Parameters across the sorted range
# This generates the "Scaling Curves"
pred_params = trained_model(z_sorted) # Shape (4, N)

# Extract individual parameters
V_max_pred = vec(pred_params[1, :])
K_pred     = vec(pred_params[2, :])
m_pred     = vec(pred_params[3, :])
Y_pred     = vec(pred_params[4, :])
z_axis     = vec(z_sorted)

# 4. Plot
# We use a 2x2 layout to see all trade-offs at once
p1 = plot(z_axis, V_max_pred, 
    title="Maximum Uptake (Vmax)", ylabel="Rate (1/h)", 
    lw=3, label=:none, color=:blue)

p2 = plot(z_axis, K_pred, 
    title="Half-Sat Constant (K)", ylabel="Concentration (mg/L)", 
    lw=3, label=:none, color=:orange)

p3 = plot(z_axis, m_pred, 
    title="Maintenance Cost (m)", ylabel="Rate (1/h)", 
    lw=3, label=:none, color=:red)

p4 = plot(z_axis, Y_pred, 
    title="Yield (Y)", ylabel="Efficiency (g/g)", 
    lw=3, label=:none, color=:green)

# Combine
final_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(800, 600), plot_title="Learned Metabolic Scaling Laws")
display(final_plot)