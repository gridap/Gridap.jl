
"""
Abstract type representing an iterable collection of Arrays{T,N},
where each array is associated to a cell.
"""
abstract type CellArray{T,N} end

Base.iterate(::CellArray)::Union{Nothing,Tuple{AbstractArray{T,N},Any}} = @abstractmethod

Base.iterate(::CellArray,state)::Union{Nothing,Tuple{AbstractArray{T,N},Any}} = @abstractmethod

Base.length(::CellArray)::Int = @abstractmethod

cellsize(::CellArray{T,N} where {T,N})::NTuple{N,Int} = @abstractmethod

Base.IteratorEltype(::Type{C} where C <: CellArray{T,N} where {T,N}) = EltypeUnknown()

cellsize(self::CellArray,i::Int) = (s = cellsize(self); s[i])

celllength(self::CellArray) = prod(cellsize(self))

function Base.show(io::IO,self::CellArray)
  for (i, a) in enumerate(self)
    println(io,"$i -> $a")
  end
end

const CellVector{T} = CellArray{T,1} where T

"""
Abstract type representing an indexable CellArray.
By implementing a concrete IndexableCellArray, one automatically
gets a type that is also iterable
"""
abstract type IndexableCellArray{T,N} <: CellArray{T,N} end

Base.getindex(::IndexableCellArray{T,N} where {T,N},cell::Int)::AbstractArray{T,N} = @abstractmethod

@inline Base.iterate(self::IndexableCellArray) = iterate(self,0)

@inline function Base.iterate(self::IndexableCellArray,state::Int)
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
abstract type CellArrayFromUnaryOp{C<:CellArray,T,N} <: CellArray{T,N} end

inputcellarray(::CellArrayFromUnaryOp{C,T,N} where {C,T,N})::C where C = @abstractmethod

computesize(::CellArrayFromUnaryOp, asize) = @abstractmethod

computevals!(::CellArrayFromUnaryOp, a, v) = @abstractmethod

Base.length(self::CellArrayFromUnaryOp) = length(inputcellarray(self))

cellsize(self::CellArrayFromUnaryOp) = computesize(self,cellsize(inputcellarray(self)))

@inline function Base.iterate(self::CellArrayFromUnaryOp{C,T,N}) where {C,T,N}
  u = Array{T,N}(undef,cellsize(self))
  v = CachedArray(u)
  anext = iterate(inputcellarray(self))
  if anext === nothing; return nothing end
  iteratekernel(self,anext,v)
end

@inline function Base.iterate(self::CellArrayFromUnaryOp,state)
  v, astate = state
  anext = iterate(inputcellarray(self),astate)
  if anext === nothing; return nothing end
  iteratekernel(self,anext,v)
end

function iteratekernel(self::CellArrayFromUnaryOp,anext,v)
  a, astate = anext
  vsize = computesize(self,size(a))
  setsize!(v,vsize)
  computevals!(self,a,v)
  state = (v, astate)
  (v,state)
end

"""
Like CellArrayFromUnaryOp but for the particular case of element-wise operation
in the elements of the returned array
"""
abstract type CellArrayFromElemUnaryOp{C,T,N} <: CellArrayFromUnaryOp{C,T,N} end

computesize(::CellArrayFromElemUnaryOp, asize) = asize

"""
Abstract type to be used for the implementation of types representing
the lazy result of applying a binary operation on two CellArray objects
"""
abstract type CellArrayFromBinaryOp{A<:CellArray,B<:CellArray,T,N} <: CellArray{T,N} end

leftcellarray(::CellArrayFromBinaryOp{A,B,T,N} where {A,B,T,N})::A = @abstractmethod

rightcellarray(::CellArrayFromBinaryOp{A,B,T,N} where {A,B,T,N})::B = @abstractmethod

computesize(::CellArrayFromBinaryOp, asize, bsize) = @abstractmethod

computevals!(::CellArrayFromBinaryOp, a, b, v) = @abstractmethod

function Base.length(self::CellArrayFromBinaryOp)
  @assert length(rightcellarray(self)) == length(leftcellarray(self))
  length(rightcellarray(self))
end

cellsize(self::CellArrayFromBinaryOp) = computesize(self,cellsize(leftcellarray(self)),cellsize(rightcellarray(self)))

@inline function Base.iterate(self::CellArrayFromBinaryOp{A,B,T,N}) where {A,B,T,N}
  u = Array{T,N}(undef,cellsize(self))
  v = CachedArray(u)
  anext = iterate(leftcellarray(self))
  if anext === nothing; return nothing end
  bnext = iterate(rightcellarray(self))
  if bnext === nothing; return nothing end
  iteratekernel(self,anext,bnext,v)
end

@inline function Base.iterate(self::CellArrayFromBinaryOp,state)
  v, astate, bstate = state
  anext = iterate(leftcellarray(self),astate)
  if anext === nothing; return nothing end
  bnext = iterate(rightcellarray(self),bstate)
  if bnext === nothing; return nothing end
  iteratekernel(self,anext,bnext,v)
end

function iteratekernel(self::CellArrayFromBinaryOp,anext,bnext,v)
  a, astate = anext
  b, bstate = bnext
  vsize = computesize(self,size(a),size(b))
  setsize!(v,vsize)
  computevals!(self,a,b,v)
  state = (v, astate, bstate)
  (v,state)
end

"""
Like CellArrayFromBinaryOp but for the particular case of element-wise operation
in the elements of the returned array
"""
abstract type CellArrayFromElemBinaryOp{A,B,T,N} <: CellArrayFromBinaryOp{A,B,T,N} end

function computesize(::CellArrayFromElemBinaryOp, asize, bsize)
  Base.Broadcast.broadcast_shape(asize,bsize)
end

