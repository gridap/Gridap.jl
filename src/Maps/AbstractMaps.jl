"""
Abstract map that takes an `M`-dim array of `S` values and returns an `N`-dim
array of `T` values
"""
abstract type Map{S,M,T,N} end

"""
Abstract field of rank `T` (e.g., scalar, vector, tensor) on a manifold of
dimension `D`
"""
const Field{D,T<:FieldValue} = Map{Point{D},1,T,1}

"""
Abstract basis for a space of fields of rank `T` (e.g., scalar, vector, tensor)
on a manifold of dimension `D`
"""
const Basis{D,T<:FieldValue} = Map{Point{D},1,T,2}

"""
Abstract geometry map
"""
const Geomap{D,Z} = Field{D,Point{Z}}

"""
Evaluate a `Map` on a set of points
"""
function evaluate!(
  this::Map{S,M,T,N},
  points::AbstractArray{S,M},
  v::AbstractArray{T,N}) where {S,M,T,N}
  @abstractmethod
end

"""
Create the gradient of a `Map`
"""
function gradient(
  this::Map{S,M,T,N})::Map{S,M,TG,N} where {S,M,T,N,TG}
  @abstractmethod
end

"""
Return dimension of the output array
"""
function return_size(
  ::Map{S,M,T,N},::NTuple{M,Int})::NTuple{N,Int} where {S,M,T,N}
  @abstractmethod
end

"""
Same as `evaluate!` but allocates output
"""
function evaluate(
  this::Map{S,M,T,N},
  points::AbstractArray{S,M}) where {S,M,T,N}
  v_size = return_size(this, size(points))
  v = Array{T,N}(undef, v_size)
  evaluate!(this,points,v)
  return v
end

