__precompile__()

module Gridap

# TODO move to TensorValues
using StaticArrays
using TensorValues
export mutable
mutable(::Type{MultiValue{S,T,N,L}}) where {S,T,N,L} = MArray{S,T,N,L}

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
