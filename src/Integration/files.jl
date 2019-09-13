
include("Quadratures.jl")
@reexport using Gridap.Quadratures

include("DuffyQuadratures.jl")
@reexport using Gridap.DuffyQuadratures

include("QuadratureFactories.jl")
@reexport using Gridap.QuadratureFactories

include("CellQuadratures.jl")
@reexport using Gridap.CellQuadratures

include("Triangulations.jl")
@reexport using Gridap.Triangulations

include("CellIntegration.jl")
@reexport using Gridap.CellIntegration

include("BoundaryDescriptors.jl")
@reexport using Gridap.BoundaryDescriptors

include("BoundaryTriangulations.jl")
@reexport using Gridap.BoundaryTriangulations

include("SkeletonTriangulations.jl")
@reexport using Gridap.SkeletonTriangulations

include("NormalVectors.jl")
@reexport using Gridap.NormalVectors

include("BoundaryCellFields.jl")
@reexport using Gridap.BoundaryCellFields

include("SkeletonCellFields.jl")
@reexport using Gridap.SkeletonCellFields
