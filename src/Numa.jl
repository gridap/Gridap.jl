__precompile__()

module Numa

using Base.Cartesian

include("Methods.jl")
include("Helpers.jl")

# CellArrays tools
include("Fields/FieldValues.jl")
include("Fields/Fields.jl")
include("Polynomials.jl")
include("Polytopes.jl")
include("Quadratures.jl")
include("CellValues/CellValues.jl")
include("CellFunctions/CellFunctions.jl")

include("Meshes.jl")
include("FESpaces/RefFEs.jl")
include("Geometry/Geometry.jl")

include("CellQuadratures.jl")
include("CellIntegration.jl")
include("Vtkio.jl")


# FESpaces tools
include("FESpaces/FESpaces.jl")

end #module Numa
