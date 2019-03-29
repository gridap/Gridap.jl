
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

"""
Assumes that outer is defined for instances of T and T as type as well
"""
function bouter(a::ConstantCellArray{T,N}, b::ConstantCellArray{S,N}) where {T,S,N}
  @assert length(a) == length(b)
  R = outer(T,S)
  s = Base.Broadcast.broadcast_shape(size(a.array),size(b.array))
  c = Array{R,N}(undef,s)
  c .= outer.(a.array,b.array)
  ConstantCellArray(c,a.length)
end

"""
Assumes that inner is defined for instances of T and T as type as well
"""
function binner(a::ConstantCellArray{T,N}, b::ConstantCellArray{T,N}) where {T,N}
  @assert length(a) == length(b)
  R = inner(T,T)
  s = Base.Broadcast.broadcast_shape(size(a.array),size(b.array))
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

function cellsum(self::ConstantCellArray{T,N};dim::Int) where {T,N}
  b = sum(self.array,dims=dim)
  s = cellsumsize(size(b),Val(dim))
  c = copy(reshape(b,s))
  ConstantCellArray(c,self.length)
end

function cellnewaxis(self::ConstantCellArray;dim::Int)
  s = [ v for v in size(self.array)]
  insert!(s,dim,1)
  shape = tuple(s...)
  c = copy(reshape(self.array,shape))
  ConstantCellArray(c,self.length)
end
