__precompile__()

module Numa

using Base.Cartesian

include("Methods.jl")
include("Helpers.jl")


include("FieldValues.jl")
include("CachedArrays.jl")
include("Geometry/Polytopes.jl")
include("Maps/Maps.jl")
include("Quadratures.jl")
include("CellValues/CellValues.jl")
include("CellMaps/CellMaps.jl")

include("Meshes.jl")
include("Polynomials.jl")
include("FESpaces/RefFEs.jl")

include("Geometry/Geometry.jl")

include("FESpaces/FESpaces.jl")

include("CellQuadratures.jl")
include("CellIntegration.jl")
include("Vtkio.jl")

end #module Numa
