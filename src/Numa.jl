__precompile__()

module Numa
using Base.Cartesian
include("Polynomial.jl")
#include("Polytope.jl")
#include("RefFE.jl")
#include("Quadrature.jl")
#include("Mesh.jl")
#include("FESpace.jl")
#include("BilinearForm.jl")
end #module
