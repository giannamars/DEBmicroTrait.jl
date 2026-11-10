loss(prediction, target, p::ModelConfig) = sum(abs2, prediction .- target)
