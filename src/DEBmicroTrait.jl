module DEBmicroTrait

using FieldMetadata,
      Unitful,
      Roots,
      LinearAlgebra,
      Statistics,
      DataFrames,
      DataFramesMeta

import FieldMetadata: @default, @description, @units, @bounds, @logscaled, @flattenable, @plottable, @selectable, @chain,
                      default, description, units, bounds, logscaled, flattenable, plottable, selectable

# Field metadata columns
@chain columns @units @description

using Unitful: J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol
Unitful.register(@__MODULE__);
@unit bp "bp" Basepair 1u"m" true;

export AbstractMetabolism, AbstractMetabolismC, MetabolismC
export AbstractAssimilation, AbstractAssimilationC, AssimilationC
export AbstractTurnover, Turnover
export AbstractParams, Params
export AbstractIsolate, IsolateComposition
export AbstractSetup, Setup

include("components/assimilation.jl")
include("components/metabolism.jl")
include("components/soil_properties.jl")
include("components/supeca.jl")
include("components/thermostoichwizard.jl")
include("components/turnover.jl")

include("allometry/cell_assimilation.jl")
include("allometry/cell_composition.jl")
include("allometry/cell_metabolism.jl")
include("allometry/cell_turnover.jl")

include("setup.jl")
include("model.jl")
end
