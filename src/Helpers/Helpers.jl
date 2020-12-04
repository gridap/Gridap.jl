"""
This module provides a set of helper macros and helper functions

The exported macros are:

$(EXPORTS)

"""
module Helpers
using DocStringExtensions

import Base: +, -, *, /, transpose, adjoint
import LinearAlgebra: cross, tr, dot

export @abstractmethod
export @notimplemented
export @notimplementedif
export @unreachable
export @check
export tfill
export get_val_parameter
export first_and_tail
export GridapType
#export operate

include("Macros.jl")

include("HelperFunctions.jl")

include("GridapTypes.jl")

end # module Helpers
