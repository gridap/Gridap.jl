
include("FESpaces.jl")
@reexport using Gridap.FESpaces

include("ConformingFESpaces.jl")
@reexport using Gridap.ConformingFESpaces

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
