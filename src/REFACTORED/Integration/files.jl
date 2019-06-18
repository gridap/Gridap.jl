
include("Quadratures.jl")
@reexport using Gridap.Quadratures

include("CellQuadratures.jl")
@reexport using Gridap.CellQuadratures

include("Triangulations.jl")
@reexport using Gridap.Triangulations

include("CellIntegration.jl")
@reexport using Gridap.CellIntegration
