__precompile__()

module Numa

using Base.Cartesian

include("Methods.jl")
include("Helpers.jl")

# CellArrays tools
include("FieldValues.jl")
include("Maps/Maps.jl")
include("Polynomials.jl")
include("Polytopes.jl")
include("Quadratures.jl")
include("CellValues/CellValues.jl")
include("CellMaps/CellMaps.jl")

include("Meshes.jl")
include("FESpaces/RefFEs.jl")
include("Geometry/Geometry.jl")

include("CellQuadratures.jl")
include("CellIntegration.jl")
include("Vtkio.jl")


# FESpaces tools
include("FESpaces/FESpaces.jl")

end #module Numa
