__precompile__()

module Gridap

using Reexport

include("Utils/files.jl")

include("CellValues/files.jl")

include("Fields/files.jl")

include("RefFEs/files.jl")

include("Integration/files.jl")

include("Geometry/files.jl")

include("Algebra/files.jl")

include("FESpaces/files.jl")

include("MultiField/files.jl")

include("Visualization/files.jl")

end # module
