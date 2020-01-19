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

include("Macros.jl")

include("HelperFunctions.jl")

include("GridapTypes.jl")

end # module Helpers
