
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

"""
    first_and_tail(a::Tuple)

Equivalent to `(first(a), Base.tail(a))`.
"""
function first_and_tail(a::Tuple)
  first(a), Base.tail(a)
end

