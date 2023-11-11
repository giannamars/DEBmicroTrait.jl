module DEBmicroTrait

using SafeTestsets,
      Roots,
      LinearAlgebra,
      Distributions,
      Dierckx

export AbstractMetabolism, AbstractMetabolismC, MetabolismC
export AbstractAssimilation, AbstractAssimilationC, AssimilationC, AssimilationCM
export AbstractTurnover, Turnover
export AbstractForcing, Forcing
export AbstractDepolymerization, Depolymerization, DepolymerizationM
export AbstractParams, Params, ParamsF
export AbstractSetup, Setup
#
include("components/assimilation.jl")
include("components/metabolism.jl")
include("components/soil_properties.jl")
include("components/supeca.jl")
include("components/thermostoichwizard.jl")
include("components/turnover.jl")
include("components/forcing.jl")
include("components/depolymerization.jl")
#
include("allometry/cell_assimilation.jl")
include("allometry/cell_composition.jl")
include("allometry/cell_metabolism.jl")
include("allometry/cell_turnover.jl")
include("allometry/cell_depolymerization.jl")
#
include("setup.jl")
include("model.jl")
end
