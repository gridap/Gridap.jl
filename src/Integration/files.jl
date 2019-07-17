
include("Quadratures.jl")
@reexport using Gridap.Quadratures

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

include("BoundaryCellFields.jl")
@reexport using Gridap.BoundaryCellFields

include("SkeletonCellFields.jl")
@reexport using Gridap.SkeletonCellFields
