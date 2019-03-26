#TODO This is a temporary file in order to explore an alternative CellArray design.

"""
Abstract type representing an iterable collection of Arrays{T,N},
where each array is associated to a cell.
"""
abstract type OtherCellArray{T,N} end

Base.iterate(::OtherCellArray)::Union{Nothing,Tuple{Array{T,N},Any}} = @abstractmethod

Base.iterate(::OtherCellArray,state)::Union{Nothing,Tuple{Array{T,N},Any}} = @abstractmethod

Base.length(::OtherCellArray)::Int = @abstractmethod

maxsize(::OtherCellArray{T,N} where {T,N})::NTuple{N,Int} = @abstractmethod

Base.eltype(::Type{C}) where C<:OtherCellArray{T,N} where {T,N} = Array{T,N}

maxsize(self::OtherCellArray,i::Int) = (s = maxsize(self); s[i])

maxlength(self::OtherCellArray) = prod(maxsize(self))

function Base.show(io::IO,self::OtherCellArray)
  for (i,a) in enumerate(self)
    println(io,"$i -> $a")
  end
end

"""
Abstract type representing an indexable CellArray.
By implementing a concrete IndexableCellArray, one automatically
gets a type that is also iterable
"""
abstract type OtherIndexableCellArray{T,N} <: OtherCellArray{T,N} end

Base.getindex(::OtherIndexableCellArray{T,N} where {T,N},cell::Int)::Array{T,N} = @abstractmethod

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

computevals!(::OtherCellArrayFromUnaryOp, a, v) = @abstractmethod

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
  a, astate = anext
  s = computesize(self,size(a))
  if size(v) != s; @notimplemented end
  computevals!(self,a,v)
  state = (v, astate)
  (v,state)
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
  self.array
end

Base.length(self::OtherConstantCellArray) = self.length

maxsize(self::OtherConstantCellArray) = size(self.array)

###struct OtherCellArrayFromUnaryOp{S,M,T,N} <: OtherCellArray{T,N}
###  a::OtherCellArray{S,M}
###  computevals!
###  computesize
###end
###
###Base.length(self::OtherCellArrayFromUnaryOp) = length(self.a)
###
###maxsize(self::OtherCellArrayFromUnaryOp) = self.computesize(maxsize(self.a))
###
###function Base.iterate(self::OtherCellArrayFromUnaryOp{S,M,T,N}) where {S,M,T,N}
###  values = Array{T,N}(undef,maxsize(self))
###  anext = iterate(self.a)
###  state = (values,anext)
###  iterate(self,state)
###end
###
###function Base.iterate(self::OtherCellArrayFromUnaryOp,state)
###  (values,anext) = state
###  if anext == nothing
###    nothing
###  else
###    ((avals,asize),astate) = anext
###    s = self.computesize(asize)
###    self.computevals!(avals,asize,values,s)
###    anext = iterate(self.a,astate)
###    state = (values,anext)
###    ((values,s), state)
###  end
###end
###
###struct DummyCellArray <: AbstractOtherCellArrayFromUnaryOp{Float64,1,Float64,2}
###  a::OtherConstantCellArray{Float64,1}
###end
###
###inputcellarray(self::DummyCellArray) = self.a
###
###computesize(self::DummyCellArray,asize) = (2,asize[1])
###
###function computevals!(self::DummyCellArray,a,as,v,vs)
###  @inbounds for i in 1:as[1]
###    v[1,i] = a[i]
###    v[2,i] = a[i]
###  end
###end

