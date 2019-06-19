# TODO to be moved to Tensor Values (begin)

using TensorValues
import Base: conj

conj(a::MultiValue) = MultiValue(conj(a.array))

# (end)

include("Utils/files.jl")

include("CellValues/files.jl")

include("Fields/files.jl")

include("RefFEs/files.jl")

include("Integration/files.jl")

include("Geometry/files.jl")

include("Algebra/files.jl")

