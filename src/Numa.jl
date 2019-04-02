__precompile__()

module Numa

using Base.Cartesian

include("Helpers.jl")

# CellArrays tools
include("FieldValues.jl")
include("Quadratures.jl")
include("Polynomials.jl")
include("NewPolynomials.jl")
include("CellArrays/CellArrays.jl")
include("CellQuadratures.jl")
include("CellFunctions/CellFunctions.jl")
include("IntegrationMeshes.jl")
include("Polytopes.jl")

# FESpaces tools
# include("FESpaces/RefFEs.jl")
# include("Meshes.jl")
# include("FESpaces/FESpaces.jl")

# include("BilinearForms.jl")

end #module Numa
