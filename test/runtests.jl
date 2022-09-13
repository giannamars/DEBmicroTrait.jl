using DEBmicroTrait
using SafeTestsets

@safetestset "DEBmicroTrait" begin
    @safetestset "Assimilation" begin
        include("test_assimilation.jl")
    end
    @safetestset "SUPECA" begin
        include("test_supeca.jl")
    end
    @safetestset "Metabolism" begin
        include("test_metabolism.jl")
    end
    @safetestset "ThermoStoichWizard" begin
        include("test_thermostoichiometry.jl")
    end
    @safetestset "Turnover" begin
        include("test_turnover.jl")
    end
    @safetestset "Setup" begin
        include("test_setup.jl")
    end
    @safetestset "Model" begin
        include("test_model.jl")
    end
end
