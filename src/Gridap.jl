__precompile__()

module Gridap

using Base.Cartesian

include("Methods.jl")
include("Helpers.jl")

include("FieldValues.jl")
include("CachedArrays.jl")
include("Maps/Maps.jl")
include("CellValues/CellValues.jl")
include("CellMaps/CellMaps.jl")

include("Geometry/Polytopes.jl")
include("Quadratures.jl")
include("Meshes.jl")
include("Polynomials.jl")
include("FESpaces/RefFEs.jl")

include("Geometry/Geometry.jl")
include("CellQuadratures.jl")
include("CellIntegration.jl")

include("Algebra/LinearSolvers.jl")
include("Algebra/NonLinearSolvers.jl")

include("FESpaces/FESpaces.jl")
include("FESpaces/Assemblers.jl")
include("FESpaces/FEOperators.jl")

include("MultiField/MultiCellArrays.jl")

include("Vtkio.jl")

end #module Gridap
