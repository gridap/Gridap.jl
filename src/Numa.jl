__precompile__()

module Numa

using Base.Cartesian

include("Helpers.jl")
include("CellArrays/CellArrays.jl")
include("Points.jl")
include("Quadratures.jl")
include("Polynomials.jl")
#include("CellArrays.jl")
#include("CellQuadratures.jl")
#include("EvaluableCellArrays.jl")
#include("CellFields.jl")
#include("CellBases.jl")
#include("CellBasesImpl.jl")
#include("CellFieldsImpl.jl")
#include("CellScalarsVectorsAndMatrices.jl")
#include("IntegrationMeshes.jl")
include("Polytopes.jl")
include("RefFEs.jl")
include("Meshes.jl")
include("FESpaces.jl")
# include("BilinearForms.jl")

end #module Numa
