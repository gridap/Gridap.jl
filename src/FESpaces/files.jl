
include("FESpaces.jl")
@reexport using Gridap.FESpaces

include("ConformingFESpaces.jl")
@reexport using Gridap.ConformingFESpaces

include("CLagrangianFESpaces.jl")
@reexport using Gridap.CLagrangianFESpaces

include("DLagrangianFESpaces.jl")
@reexport using Gridap.DLagrangianFESpaces

include("DiscFESpaces.jl")
@reexport using Gridap.DiscFESpaces

include("ConstrainedFESpaces.jl")
@reexport using Gridap.ConstrainedFESpaces

include("FEFunctions.jl")
@reexport using Gridap.FEFunctions

include("FEBases.jl")
@reexport using Gridap.FEBases

include("Assemblers.jl")
@reexport using Gridap.Assemblers

include("FETerms.jl")
@reexport using Gridap.FETerms

include("FEOperators.jl")
@reexport using Gridap.FEOperators
