__precompile__()

module Numa

using Base.Cartesian

include("Helpers.jl")

# CellArrays tools
include("FieldValues.jl")
include("Fields.jl")
include("Fields2.jl")#@fverdugo alternative implementation proposal
include("Quadratures.jl")
include("Polynomials.jl")
# include("NewPolynomials.jl")
include("CellValues/CellValues.jl")
include("CellFunctions/CellFunctions.jl")
include("CellQuadratures.jl")
include("CellIntegration.jl")
include("Polytopes.jl")

# FESpaces tools
include("FESpaces/RefFEs.jl")
# include("Meshes.jl")
# include("FESpaces/FESpaces.jl")

# include("BilinearForms.jl")

end #module Numa
