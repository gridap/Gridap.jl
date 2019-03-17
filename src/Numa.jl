__precompile__()

module Numa
using Base.Cartesian
include("Helpers.jl")
include("CellArrays.jl")
include("CellQuadratures.jl")
include("Polynomials.jl")
include("Quadratures.jl")
#include("Polytope.jl")
#include("RefFE.jl")
#include("Mesh.jl")
#include("FESpace.jl")
#include("BilinearForm.jl")
end #module
