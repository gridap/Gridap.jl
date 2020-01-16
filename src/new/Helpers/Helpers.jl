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
export GridapType
export get_val_parameter

"""
    @abstractmethod

Macro used in generic functions that must be overloaded by derived types.
"""
macro abstractmethod()
  quote
    error("This function belongs to an interface definition and cannot be used.")
  end
end

"""
    @notimplemented
    @notimplemented "Error message"

Macro used to raise an error, when something is not implemented.
"""
macro notimplemented(message="This function is not yet implemented")
  quote
    error($(esc(message)))
  end
end

"""
    @notimplementedif condition
    @notimplementedif condition "Error message"

Macro used to raise an error if the `condition` is true
"""
macro notimplementedif(condition,message="This function is not yet implemented")
  quote
    if $(esc(condition))
      @notimplemented $(esc(message))
    end
  end
end

"""
    @unreachable
    @unreachable "Error message"

Macro used to make sure that a line of code is never reached.
"""
macro unreachable(message="This line of code cannot be reached")
  quote
    error($(esc(message)))
  end
end

"""
    tfill(v, ::Val{D}) where D

Returns a tuple of length `D` that contains `D` times the object `v`.
In contrast to `tuple(fill(v,D)...)` which returns the same result, this function is type-stable.
"""
function tfill(v, ::Val{D}) where D
  t = tfill(v, Val{D-1}())
  (v,t...)
end

tfill(v,::Val{0}) = ()
tfill(v,::Val{1}) = (v,)
tfill(v,::Val{2}) = (v,v)
tfill(v,::Val{3}) = (v,v,v)


"""
    abstract type GridapType end
"""
abstract type GridapType end

function show(io::IO,object::GridapType)
  print(io,"$(nameof(typeof(object)))()")
end


"""
    get_val_parameter(::Val{T}) where T
    get_val_parameter(::Type{Val{T}}) where T

Returns `T`.
"""
function get_val_parameter(::Val{T}) where T
  T
end

function get_val_parameter(::Type{Val{T}}) where T
  T
end

end # module Helpers
