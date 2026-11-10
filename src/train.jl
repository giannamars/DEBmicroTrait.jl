function training_step!(model::SoilModel, batch; mass_penalty = 0.0)
    prediction, target = batch
    return loss(prediction, target, model; mass_penalty = mass_penalty)
end

