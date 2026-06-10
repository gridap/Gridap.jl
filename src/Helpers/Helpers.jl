"""
This module provides a set of helper macros and helper functions

$(public_names_in_md(@__MODULE__))
"""
module Helpers
using DocStringExtensions
using Test

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
export @check_inferred
export tfill
export get_val_parameter
export first_and_tail
export GridapType
export GridapLocalInt, set_local_integer_type
export set_debug_mode, set_performance_mode
export default_num_nearest_vertices, set_num_nearest_vertices
export public_names_in_md

include("Preferences.jl")

include("Macros.jl")

include("HelperFunctions.jl")

include("GridapTypes.jl")

end # module Helpers
