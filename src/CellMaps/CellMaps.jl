module CellMaps

using Reexport

include("AbstractCellMaps.jl")

@reexport using Gridap.CellMaps.AbstractCellMaps

include("CellMapValues.jl")

include("Operators.jl")

include("Composition.jl")

include("CellBasesWithGeomap.jl")

include("ConstantCellMaps.jl")

include("Testers.jl")

end #module CellMaps

