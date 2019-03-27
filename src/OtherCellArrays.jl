module OtherCellArrays

using LinearAlgebra: det
import LinearAlgebra

using Numa.Helpers

export OtherCellArray
export IndexableCellArray
export OtherCellArrayFromUnaryOp
export OtherCellArrayFromElemUnaryOp
export OtherConstantCellArray
export maxsize
export maxlength

#TODO This is a temporary file in order to explore an alternative CellArray design.

"""
Abstract type representing an iterable collection of Arrays{T,N},
where each array is associated to a cell.
"""
abstract type OtherCellArray{T,N} end

Base.iterate(::OtherCellArray)::Union{Nothing,Tuple{Tuple{Array{T,N},NTuple{N,Int}},Any}} = @abstractmethod

Base.iterate(::OtherCellArray,state)::Union{Nothing,Tuple{Tuple{Array{T,N},NTuple{N,Int}},Any}} = @abstractmethod

Base.length(::OtherCellArray)::Int = @abstractmethod

maxsize(::OtherCellArray{T,N} where {T,N})::NTuple{N,Int} = @abstractmethod

"""
Abstract type representing an indexable CellArray.
By implementing a concrete IndexableCellArray, one automatically
gets a type that is also iterable
"""
abstract type OtherIndexableCellArray{T,N} <: OtherCellArray{T,N} end

Base.getindex(::OtherIndexableCellArray{T,N} where {T,N},cell::Int)::Tuple{Array{T,N},NTuple{N,Int}} = @abstractmethod

"""
Abstract type to be used for the implementation of types representing
the lazy result of applying an unary operation on a CellArray
"""
abstract type OtherCellArrayFromUnaryOp{C<:OtherCellArray,T,N} <: OtherCellArray{T,N} end

inputcellarray(::OtherCellArrayFromUnaryOp{C,T,N} where {C,T,N})::C = @abstractmethod

computesize(::OtherCellArrayFromUnaryOp, asize) = @abstractmethod

computevals!(::OtherCellArrayFromUnaryOp, a, asize, v, vsize) = @abstractmethod

"""
Like OtherCellArrayFromUnaryOp but for the particular case of element-wise operation
in the elements of the returned array
"""
abstract type OtherCellArrayFromElemUnaryOp{C,T,N} <: OtherCellArrayFromUnaryOp{C,T,N} end

# Concrete implementations

"""
Concrete implementation of CellArray, where the same array
is associated to all cells. Typically, this is useful for
discretizations with a single cell type.
"""
struct OtherConstantCellArray{T,N} <: OtherIndexableCellArray{T,N}
  array::Array{T,N}
  length::Int
end

"""
Type that stores the lazy result of evaluating the determinant
of each element in a CellArray
"""
struct OtherCellArrayFromDet{C,T,N} <: OtherCellArrayFromElemUnaryOp{C,T,N}
  a::C
end

"""
Type that stores the lazy result of evaluating the inverse of
of each element in a CellArray
"""
struct OtherCellArrayFromInv{C,T,N} <: OtherCellArrayFromElemUnaryOp{C,T,N}
  a::C
end

# Methods

# OtherCellArray

Base.eltype(::Type{C}) where C<:OtherCellArray{T,N} where {T,N} = Array{T,N}

maxsize(self::OtherCellArray,i::Int) = (s = maxsize(self); s[i])

maxlength(self::OtherCellArray) = prod(maxsize(self))

function Base.show(io::IO,self::OtherCellArray)
  for (i,(a,s)) in enumerate(self)
    v = viewtosize(a,s)
    println(io,"$i -> $v")
  end
end

function Base.:(==)(a::OtherCellArray{T,N},b::OtherCellArray{T,N}) where {T,N}
  length(a) != length(b) && return false
  maxsize(a) != maxsize(b) && return false
  if N != 1; @notimplemented end
  for ((ai,ais),(bi,bis)) in zip(a,b)
    ais != bis && return false
    for j in 1:ais[1]
      ai[j] != bi[j] && return false
    end
  end
  return true
end

"""
Assumes that det is defined for instances of T
and that the result is Float64
"""
function LinearAlgebra.det(self::OtherCellArray{T,N}) where {T,N}
  OtherCellArrayFromDet{typeof(self),Float64,N}(self)
end

"""
Assumes that inv is defined for instances of T
"""
function LinearAlgebra.inv(self::OtherCellArray{T,N}) where {T,N}
  OtherCellArrayFromInv{typeof(self),T,N}(self)
end

# OtherIndexableCellArray

Base.iterate(self::OtherIndexableCellArray) = iterate(self,0)

function Base.iterate(self::OtherIndexableCellArray,state::Int)
  if length(self) == state
    nothing
  else
    k = state+1
    (self[k],k)
  end
end

# OtherCellArrayFromUnaryOp

Base.length(self::OtherCellArrayFromUnaryOp) = length(inputcellarray(self))

maxsize(self::OtherCellArrayFromUnaryOp) = computesize(self,maxsize(inputcellarray(self)))

@inline function Base.iterate(self::OtherCellArrayFromUnaryOp{C,T,N}) where {C,T,N}
  v = Array{T,N}(undef,maxsize(self))
  anext = iterate(inputcellarray(self))
  if anext === nothing; return nothing end
  iteratekernel(self,anext,v)
end

@inline function Base.iterate(self::OtherCellArrayFromUnaryOp,state)
  v, astate = state
  anext = iterate(inputcellarray(self),astate)
  if anext === nothing; return nothing end
  iteratekernel(self,anext,v)
end

function iteratekernel(self::OtherCellArrayFromUnaryOp,anext,v)
  (a,asize), astate = anext
  vsize = computesize(self,asize)
  computevals!(self,a,asize,v,vsize)
  state = (v, astate)
  ((v,vsize),state)
end

# OtherCellArrayFromElemUnaryOp

computesize(::OtherCellArrayFromElemUnaryOp, asize) = asize

# OtherConstantCellArray

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

# OtherCellArrayFromDet

inputcellarray(self::OtherCellArrayFromDet) = self.a

function computevals!(::OtherCellArrayFromDet, a, asize, v, vsize)
  if length(asize) != 1; @notimplemented end
  for i in 1:asize[1]
    v[i] = det(a[i])
  end
end

# OtherCellArrayFromInv

inputcellarray(self::OtherCellArrayFromInv) = self.a

function computevals!(::OtherCellArrayFromInv, a, asize, v, vsize)
  if length(asize) != 1; @notimplemented end
  for i in 1:asize[1]
    v[i] = inv(a[i])
  end
end

end # module OtherCellArrays
