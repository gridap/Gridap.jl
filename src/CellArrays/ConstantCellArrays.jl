
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

function Base.:+(a::ConstantCellArray{T,N},b::ConstantCellArray{T,N}) where {T,N}
  @assert size(a.array) == size(b.array)
  @assert length(a) == length(b)
  c = Array{T,N}(undef,size(a.array))
  c .= a.array .+ b.array
  ConstantCellArray(c,a.length)
end

function Base.:-(a::ConstantCellArray{T,N},b::ConstantCellArray{T,N}) where {T,N}
  @assert size(a.array) == size(b.array)
  @assert length(a) == length(b)
  c = Array{T,N}(undef,size(a.array))
  c .= a.array .- b.array
  ConstantCellArray(c,a.length)
end

function Base.:*(a::ConstantCellArray{T,N},b::ConstantCellArray{T,N}) where {T,N}
  @assert size(a.array) == size(b.array)
  @assert length(a) == length(b)
  c = Array{T,N}(undef,size(a.array))
  c .= a.array .* b.array
  ConstantCellArray(c,a.length)
end

function Base.:/(a::ConstantCellArray{T,N},b::ConstantCellArray{T,N}) where {T,N}
  @assert size(a.array) == size(b.array)
  @assert length(a) == length(b)
  c = Array{T,N}(undef,size(a.array))
  c .= a.array ./ b.array
  ConstantCellArray(c,a.length)
end

"""
Assumes that det is defined for instances of T
and that the result is Float64
"""
function LinearAlgebra.det(self::ConstantCellArray{T,N}) where {T,N}
  deta = Array{Float64,N}(undef,size(self.array))
  deta .= det.(self.array)
  ConstantCellArray(deta,self.length)
end

"""
Assumes that inv is defined for instances of T
"""
function LinearAlgebra.inv(self::ConstantCellArray{T,N}) where {T,N}
  deta = Array{T,N}(undef,size(self.array))
  deta .= inv.(self.array)
  ConstantCellArray(deta,self.length)
end
