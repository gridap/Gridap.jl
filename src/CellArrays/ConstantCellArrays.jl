
"""
Concrete implementation of CellArray, where the same array
is associated to all cells. Typically, this is useful for
discretizations with a single cell type.
"""
struct ConstantCellArray{T,N} <: IndexableCellArray{T,N}
  array::Array{T,N}
  length::Int
end

function Base.getindex(self::ConstantCellArray,cell::Int)
  @assert 1 <= cell
  @assert cell <= length(self)
  self.array
end

Base.length(self::ConstantCellArray) = self.length

Base.eltype(::Type{ConstantCellArray{T,N}}) where {T,N} = Array{T,N}

cellsize(self::ConstantCellArray) = size(self.array)

function Base.:(==)(a::ConstantCellArray{T,N},b::ConstantCellArray{T,N}) where {T,N}
  a.array != b.array && return false
  a.length != b.length && return false
  return true
end

function Base.:+(a::ConstantCellArray,b::ConstantCellArray)
  @assert length(a) == length(b)
  c = a.array .+ b.array
  ConstantCellArray(c,a.length)
end

function Base.:-(a::ConstantCellArray,b::ConstantCellArray)
  @assert length(a) == length(b)
  c = a.array .- b.array
  ConstantCellArray(c,a.length)
end

function Base.:*(a::ConstantCellArray,b::ConstantCellArray)
  @assert length(a) == length(b)
  c = a.array .* b.array
  ConstantCellArray(c,a.length)
end

function Base.:/(a::ConstantCellArray,b::ConstantCellArray)
  @assert length(a) == length(b)
  c = a.array ./ b.array
  ConstantCellArray(c,a.length)
end

function outer(a::ConstantCellArray{T,N} where T, b::ConstantCellArray{S,N} where S) where N
  @assert length(a) == length(b)
  R = outer(T,S)
  s = Base.Broadcast.broadcast_shape(size(a),size(b))
  c = Array{R,N}(undef,s)
  c .= outer.(a.array,b.array)
  ConstantCellArray(c,a.length)
end

function inner(a::ConstantCellArray{T,N}, b::ConstantCellArray{T,N}) where {T,N}
  @assert length(a) == length(b)
  R = inner(T,S)
  s = Base.Broadcast.broadcast_shape(size(a),size(b))
  c = Array{R,N}(undef,s)
  c .= inner.(a.array,b.array)
  ConstantCellArray(c,a.length)
end

"""
Assumes that det is defined for instances of T
and that the result is Float64
"""
function LinearAlgebra.det(self::ConstantCellArray)
  deta = det.(self.array)
  ConstantCellArray(deta,self.length)
end

"""
Assumes that inv is defined for instances of T
"""
function LinearAlgebra.inv(self::ConstantCellArray)
  deta = inv.(self.array)
  ConstantCellArray(deta,self.length)
end

function cellsum(self::ConstantCellArray{T,N};dims::Int) where {T,N}
  b = sum(self.array,dims=dims)
  @notimplementedif dims != N
  sb = size(b)
  s = tuple([v for (i,v) in enumerate(sb) if i<length(sb) ]...)
  c = copy(reshape(b,s))
  ConstantCellArray(c,self.length)
end

function cellreshape(self::ConstantCellArray,shape::NTuple{M,Int}) where M
  c = copy(reshape(self.array,shape))
  ConstantCellArray(c,self.length)
end
