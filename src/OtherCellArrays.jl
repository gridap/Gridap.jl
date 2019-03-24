#TODO This is a temporary file in order to explore an alternative CellArray design.
# We need to carefully think and decide which version in better since this is a very
# core component, which will affect all the project.

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
    #TODO: We need to restrict a to s
    println(io,"$i -> $a")
  end
end

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

abstract type AbstractOtherCellArrayFromUnaryOp{S,M,T,N} <: OtherCellArray{T,N} end

inputcellarray(
  ::AbstractOtherCellArrayFromUnaryOp{S,M,T,N} where {S,M,T,N})::CellArray{S,M} = @abstractmethod

computesize(::AbstractOtherCellArrayFromUnaryOp, asize) = @abstractmethod

computevals!(::AbstractOtherCellArrayFromUnaryOp, a, asize ,v, vsize) = @abstractmethod

Base.length(self::AbstractOtherCellArrayFromUnaryOp) = length(inputcellarray(self))

maxsize(self::AbstractOtherCellArrayFromUnaryOp) = computesize(self,maxsize(inputcellarray(self)))

function Base.iterate(self::AbstractOtherCellArrayFromUnaryOp{S,M,T,N}) where {S,M,T,N}
  values = Array{T,N}(undef,maxsize(self))
  anext = iterate(inputcellarray(self))
  state = (values,anext)
  iterate(self,state)
end

function Base.iterate(self::AbstractOtherCellArrayFromUnaryOp,state)
  (values,anext) = state
  if anext == nothing
    nothing
  else
    ((avals,asize),astate) = anext
    s = computesize(self,asize)
    computevals!(self,avals,asize,values,s)
    anext = iterate(inputcellarray(self),astate)
    state = (values,anext)
    ((values,s), state)
  end
end

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

struct OtherCellArrayFromUnaryOp{S,M,T,N} <: OtherCellArray{T,N}
  a::OtherCellArray{S,M}
  computevals!
  computesize
end

Base.length(self::OtherCellArrayFromUnaryOp) = length(self.a)

maxsize(self::OtherCellArrayFromUnaryOp) = self.computesize(maxsize(self.a))

function Base.iterate(self::OtherCellArrayFromUnaryOp{S,M,T,N}) where {S,M,T,N}
  values = Array{T,N}(undef,maxsize(self))
  anext = iterate(self.a)
  state = (values,anext)
  iterate(self,state)
end

function Base.iterate(self::OtherCellArrayFromUnaryOp,state)
  (values,anext) = state
  if anext == nothing
    nothing
  else
    ((avals,asize),astate) = anext
    s = self.computesize(asize)
    self.computevals!(avals,asize,values,s)
    anext = iterate(self.a,astate)
    state = (values,anext)
    ((values,s), state)
  end
end

struct DummyCellArray <: AbstractOtherCellArrayFromUnaryOp{Float64,1,Float64,2}
  a::OtherConstantCellArray{Float64,1}
end

inputcellarray(self::DummyCellArray) = self.a

computesize(self::DummyCellArray,asize) = (2,asize[1])

function computevals!(self::DummyCellArray,a,as,v,vs)
  @inbounds for i in 1:as[1]
    v[1,i] = a[i]
    v[2,i] = a[i]
  end
end

