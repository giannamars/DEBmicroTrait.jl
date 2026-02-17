using Flux
using Statistics, Random
using CSV, DataFrames
using MultivariateStats, DelimitedFiles
using Plots

# ==============================================================================
# 1. DATA LOADING & PREPROCESSING
# ==============================================================================

include("mags_2traits.jl") # This script loads and preprocesses the data, resulting in X_continuous and X_binary matrices ready for Flux.

# ==============================================================================
# 2. MODEL ARCHITECTURE
# ==============================================================================

struct MultiModalAE
    encoder::Chain
    shared_decoder::Chain
    head_continuous::Dense
    head_binary::Dense
end

# Register layer with Flux
Flux.@layer MultiModalAE

# Constructor: The "Nano" Architecture for Small Data (N=120, P=25)
function MultiModalAE(n_continuous::Int, n_binary::Int, latent_dim::Int)
    input_dim = n_continuous + n_binary
    
    # Hidden Layer Size: Rule of thumb is sqrt(input_dim * output_dim)
    # sqrt(25 * 2) ≈ 7 or 8 neurons.
    hidden_dim = 8 
    
    return MultiModalAE(
        # Encoder: 25 -> 8 -> 2
        Chain(
            Dense(input_dim => hidden_dim, relu),
            # High Dropout because the data is small
            Dropout(0.2),                 
            Dense(hidden_dim => latent_dim, relu)
        ),
        # Shared Decoder: 2 -> 8
        Chain(
            Dense(latent_dim => hidden_dim, relu)
        ),
        # Heads: 8 -> Outputs
        Dense(hidden_dim => n_continuous), 
        Dense(hidden_dim => n_binary)
    )
end

function (m::MultiModalAE)(x_continuous, x_binary)
    x_combined = vcat(x_continuous, x_binary)
    z = m.encoder(x_combined)
    h = m.shared_decoder(z)
    return m.head_continuous(h), m.head_binary(h)
end

# ==============================================================================
# 3. INITIALIZATION & OPTIMIZER
# ==============================================================================

# Dimensions
n_cont_feats = size(X_continuous, 1)
n_bin_feats  = size(X_binary, 1)
latent_size  = 2 

# 1. Instantiate Model
model = MultiModalAE(n_cont_feats, n_bin_feats, latent_size)

# 2. Setup Optimizer (Weight Decay + Adam)
learning_rate = 0.002
l2_lambda     = 0.05  # strong regularization

opt_state = Flux.setup(
    Flux.OptimiserChain(
        Flux.WeightDecay(l2_lambda), # Handles L2 penalty automatically
        Flux.Adam(learning_rate)
    ), 
    model
)

# 3. Loss Function (Clean: No manual L2 loop needed)
function compute_loss(model, x_c, x_b)
    pred_c, pred_b_logits = model(x_c, x_b)
    
    loss_mse = Flux.mse(pred_c, x_c)
    loss_bce = Flux.logitbinarycrossentropy(pred_b_logits, x_b)
    
    return loss_mse + loss_bce 
end

println("Model & Optimizer Initialized.")

# ==============================================================================
# 4. TRAINING LOOP
# ==============================================================================

# Data Split
n_samples = size(X_continuous, 2)
n_train   = floor(Int, 0.8 * n_samples)
indices   = shuffle(1:n_samples)

train_idxs = indices[1:n_train]
test_idxs  = indices[n_train+1:end]

x_c_train = X_continuous[:, train_idxs]
x_b_train = X_binary[:, train_idxs]
x_c_test  = X_continuous[:, test_idxs]
x_b_test  = X_binary[:, test_idxs]

println("Training on $(length(train_idxs)) samples. Testing on $(length(test_idxs)).")

epochs = 1000

train_loss_history = Float64[]
test_loss_history  = Float64[]
test_mse_history   = Float64[]
test_bce_history   = Float64[]
epoch_indices      = Int[]

for epoch in 1:epochs
    # --- A. Optimization Step (TRAIN DATA ONLY) ---
    grads = Flux.gradient(model) do m
        # 1. Generate Gaussian Noise (Float32 to match model)
        # 0.2f0 is the standard deviation (Magnitude of noise)
        noise = 0.2f0 .* randn(Float32, size(x_c_train))
        
        # 2. Add Noise to the Continuous Input
        # The model sees "corrupted" data...
        x_c_noisy = x_c_train .+ noise
        
        # 3. Predict using the NOISY input
        pred_c, pred_b_logits = m(x_c_noisy, x_b_train)
        
        # 4. Calculate Loss against the CLEAN Target
        # ...but must learn to predict the clean original data.
        loss_mse = Flux.mse(pred_c, x_c_train) 
        loss_bce = Flux.logitbinarycrossentropy(pred_b_logits, x_b_train)
        
        return loss_mse + loss_bce
    end
    
    # Update weights
    Flux.update!(opt_state, model, grads[1])
    
    # --- B. Evaluation Step ---
    if epoch % 100 == 0 
        # 1. Turn off Dropout for accurate evaluation
        Flux.testmode!(model) 
        
        # 2. Train Metrics (Clean Data, No Noise)
        p_c_tr, p_b_tr = model(x_c_train, x_b_train)
        # Note: We evaluate against clean data to see actual fit
        loss_tr = Flux.mse(p_c_tr, x_c_train) + Flux.logitbinarycrossentropy(p_b_tr, x_b_train)
        
        # 3. Test Metrics (Clean Data, New Genomes)
        p_c_te, p_b_te = model(x_c_test, x_b_test)
        mse_te = Flux.mse(p_c_te, x_c_test)
        bce_te = Flux.logitbinarycrossentropy(p_b_te, x_b_test)
        loss_te = mse_te + bce_te
        
        # 4. Turn Dropout back ON for the next training loop
        Flux.trainmode!(model) 

        push!(epoch_indices, epoch)
        push!(train_loss_history, loss_tr)
        push!(test_loss_history, loss_te)
        push!(test_mse_history, mse_te)
        push!(test_bce_history, bce_te)
        
        println("Epoch $epoch:")
        println("  Train: $(round(loss_tr, digits=4))")
        println("  Test:  $(round(loss_te, digits=4)) (MSE=$(round(mse_te, digits=4)) | BCE=$(round(bce_te, digits=4)))")
        
        if loss_te > (loss_tr * 1.2)
             println("  [WARNING] Overfitting detected.")
        end
    end
end

p1 = plot(epoch_indices, train_loss_history, label="Train Loss", lw=2, title="Total Loss")
plot!(p1, epoch_indices, test_loss_history, label="Test Loss", lw=2, linestyle=:dash)
plot!(p1, epoch_indices, test_mse_history, label="Test MSE", lw=2, linestyle=:dot)
plot!(p1, epoch_indices, test_bce_history, label="Test BCE", lw=2, linestyle=:solid)

# ==============================================================================
# 5. LATENT SPACE VISUALIZATION
# ==============================================================================
println("--- Visualizing Latent Space ---")

# 1. Turn off Dropout/Noise for extraction (Critical!)
Flux.testmode!(model) 

# 2. Extract Z for the FULL dataset
# We recombine the full dataset since we are done with the Train/Test split
full_continuous = Float32.(X_continuous) 
full_binary     = Float32.(X_binary)
full_input      = vcat(full_continuous, full_binary)

# Pass through encoder
Z_total = model.encoder(full_input) # Shape: (latent_dim, n_samples)

println("Latent Space Extracted. Shape: ", size(Z_total))

# 3. Save for DEB Model
writedlm("latent_space_Z_nano_tm1.csv", Z_total', ',') 
println("Saved 'latent_space_Z_nano_tm1.csv'")

# 4. Plotting
# Since latent_dim is already 2 (from the Nano model), we don't even need PCA!
# We can plot Z[1,:] vs Z[2,:] directly.

# We define x and y coordinates
z_x = Z_total[1, :]
z_y = Z_total[2, :]

# Plot colored by Genome Size
p = scatter(
    z_x, 
    z_y, 
    zcolor = genome_sizes, # This variable should still be in memory
    title = "Nano-Autoencoder Latent Space",
    label = "Genomes",
    xlabel = "Latent Dimension 1",
    ylabel = "Latent Dimension 2",
    c = :viridis,
    markerstrokewidth = 0,
    alpha = 0.8,
    size = (600, 500)
)

display(p)

# 1. Take your Nano-Model Z (2 x N)
# Z_total from the previous step

# 2. PCA Rotation to 1Dß
pca_1d = fit(PCA, Z_total; maxoutdim=1)
Z_score = MultivariateStats.transform(pca_1d, Z_total) # Shape (1, N)

# 3. Check direction
# PCA direction is arbitrary. We want High Z = Large Genome.
# Calculate correlation. If negative, flip the sign of Z.
correlation = cor(vec(Z_score), genome_sizes)
if correlation < 0
    Z_score = -Z_score
    println("Flipped Z-score to align with Genome Size.")
end

println("Final Complexity Score extracted.")
println("Correlation with Genome Size: ", abs(correlation))

# 4. Save
writedlm("genome_complexity_score_tm1.csv", Z_score', ',')