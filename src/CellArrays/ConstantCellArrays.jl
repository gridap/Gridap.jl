
"""
Concrete implementation of CellArray, where the same array
is associated to all cells. Typically, this is useful for
discretizations with a single cell type.
"""
struct OtherConstantCellArray{T,N} <: OtherIndexableCellArray{T,N}
  array::Array{T,N}
  length::Int
end

function Base.getindex(self::OtherConstantCellArray,cell::Int)
  @assert 1 <= cell
  @assert cell <= length(self)
  (self.array, size(self.array))
end

Base.length(self::OtherConstantCellArray) = self.length

maxsize(self::OtherConstantCellArray) = size(self.array)

function Base.:(==)(a::OtherConstantCellArray{T,N},b::OtherConstantCellArray{T,N}) where {T,N}
  a.array != b.array && return false
  a.length != b.length && return false
  return true
end

function Base.:+(a::OtherConstantCellArray{T,N},b::OtherConstantCellArray{T,N}) where {T,N}
  @assert size(a.array) == size(b.array)
  @assert length(a) == length(b)
  c = Array{T,N}(undef,size(a.array))
  c .= a.array .+ b.array
  OtherConstantCellArray(c,a.length)
end

function Base.:-(a::OtherConstantCellArray{T,N},b::OtherConstantCellArray{T,N}) where {T,N}
  @assert size(a.array) == size(b.array)
  @assert length(a) == length(b)
  c = Array{T,N}(undef,size(a.array))
  c .= a.array .- b.array
  OtherConstantCellArray(c,a.length)
end

function Base.:*(a::OtherConstantCellArray{T,N},b::OtherConstantCellArray{T,N}) where {T,N}
  @assert size(a.array) == size(b.array)
  @assert length(a) == length(b)
  c = Array{T,N}(undef,size(a.array))
  c .= a.array .* b.array
  OtherConstantCellArray(c,a.length)
end

function Base.:/(a::OtherConstantCellArray{T,N},b::OtherConstantCellArray{T,N}) where {T,N}
  @assert size(a.array) == size(b.array)
  @assert length(a) == length(b)
  c = Array{T,N}(undef,size(a.array))
  c .= a.array ./ b.array
  OtherConstantCellArray(c,a.length)
end

"""
Assumes that det is defined for instances of T
and that the result is Float64
"""
function LinearAlgebra.det(self::OtherConstantCellArray{T,N}) where {T,N}
  deta = Array{Float64,N}(undef,size(self.array))
  deta .= det.(self.array)
  OtherConstantCellArray(deta,self.length)
end

"""
Assumes that inv is defined for instances of T
"""
function LinearAlgebra.inv(self::OtherConstantCellArray{T,N}) where {T,N}
  deta = Array{T,N}(undef,size(self.array))
  deta .= inv.(self.array)
  OtherConstantCellArray(deta,self.length)
end
