module CellValues

using Numa.Helpers

import Base: iterate
import Base: length
import Base: eltype
import Base: size
import Base: getindex
import Base: IndexStyle
import Base: +, -, *, /
import LinearAlgebra: inv, det

include("Interfaces.jl")

end # module CellValues
