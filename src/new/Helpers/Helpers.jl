"""
This module provides a set of helper macros and helper functions

The exported macros are:

$(EXPORTS)

"""
module Helpers
using DocStringExtensions

export @abstractmethod
export @notimplemented
export @notimplementedif
export @unreachable
export tfill
export get_val_parameter
export GridapType
export unary_operation
export get_binary_operation
export convert_to_operable

import Base: +, -, *, transpose, adjoint
import LinearAlgebra: cross, tr, dot

include("Macros.jl")

include("HelperFunctions.jl")

include("GridapTypes.jl")

end # module Helpers
