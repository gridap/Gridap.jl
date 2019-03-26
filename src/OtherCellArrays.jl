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

Base.eltype(::Type{C}) where C<:OtherCellArray{T,N} where {T,N} = Array{T,N}

maxsize(self::OtherCellArray,i::Int) = (s = maxsize(self); s[i])

maxlength(self::OtherCellArray) = prod(maxsize(self))

function Base.show(io::IO,self::OtherCellArray)
  for (i,(a,s)) in enumerate(self)
    v = viewtosize(a,s)
    println(io,"$i -> $v")
  end
end

"""
Abstract type representing an indexable CellArray.
By implementing a concrete IndexableCellArray, one automatically
gets a type that is also iterable
"""
abstract type OtherIndexableCellArray{T,N} <: OtherCellArray{T,N} end

Base.getindex(::OtherIndexableCellArray{T,N} where {T,N},cell::Int)::Tuple{Array{T,N},NTuple{N,Int}} = @abstractmethod

Base.iterate(self::OtherIndexableCellArray) = iterate(self,0)

function Base.iterate(self::OtherIndexableCellArray,state::Int)
  if length(self) == state
    nothing
  else
    k = state+1
    (self[k],k)
  end
end

"""
Abstract type to be used for the implementation of types representing
the lazy result of applying an unary operation on a CellArray
"""
abstract type OtherCellArrayFromUnaryOp{C<:OtherCellArray,T,N} <: OtherCellArray{T,N} end

inputcellarray(::OtherCellArrayFromUnaryOp{C,T,N} where {C,T,N})::C = @abstractmethod

computesize(::OtherCellArrayFromUnaryOp, asize) = @abstractmethod

computevals!(::OtherCellArrayFromUnaryOp, a, asize, v, vsize) = @abstractmethod

Base.length(self::OtherCellArrayFromUnaryOp) = length(inputcellarray(self))

maxsize(self::OtherCellArrayFromUnaryOp) = computesize(self,maxsize(inputcellarray(self)))

function Base.iterate(self::OtherCellArrayFromUnaryOp{C,T,N}) where {C,T,N}
  v = Array{T,N}(undef,maxsize(self))
  anext = iterate(inputcellarray(self))
  if anext === nothing; return nothing end
  iteratekernel(self,anext,v)
end

function Base.iterate(self::OtherCellArrayFromUnaryOp,state)
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

"""
Like OtherCellArrayFromUnaryOp but for the particular case of element-wise operation
in the elements of the returned array
"""
abstract type OtherCellArrayFromElemUnaryOp{C,T,N} <: OtherCellArrayFromUnaryOp{C,T,N} end

computesize(::OtherCellArrayFromElemUnaryOp, asize) = asize

# Concrete implementations

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

