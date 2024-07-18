"""
This module provides a set of helper macros and helper functions

The exported macros are:

$(EXPORTS)

"""
module Helpers
using DocStringExtensions

@static if VERSION >= v"1.6"
  using Preferences
end

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
export set_debug_mode, set_performance_mode
#export operate

include("Preferences.jl")

include("Macros.jl")

include("HelperFunctions.jl")

include("GridapTypes.jl")

end # module Helpers
