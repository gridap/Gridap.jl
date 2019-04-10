__precompile__()

module Numa

using Base.Cartesian

include("Methods.jl")
include("Helpers.jl")

# CellArrays tools
include("Fields/FieldValues.jl")
include("Fields/Fields.jl")
include("Quadratures.jl")
include("Polynomials.jl")
include("CellValues/CellValues.jl")
include("CellFunctions/CellFunctions.jl")
include("CellQuadratures.jl")
include("CellIntegration.jl")
include("Polytopes.jl")

include("Geometry/Geometry.jl")

# FESpaces tools
include("FESpaces/RefFEs.jl")
# include("Meshes.jl")
# include("FESpaces/FESpaces.jl")

# include("BilinearForms.jl")

end #module Numa
