export CellArray, IndexableCellArray
export maxlength, maxsize, cellsum
export CellFieldValues, CellBasisValues
export CellValues, CellPoints
export CellScalars, CellVectors, CellMatrices
export ConstantCellArray

"""
Abstract type representing an iterable collection of Arrays{T,N},
where each array is associated to a cell.
"""
abstract type CellArray{T,N} end

Base.iterate(::CellArray)::Union{Nothing,Tuple{Array{T,N},Any}} = @abstractmethod

Base.iterate(::CellArray,state)::Union{Nothing,Tuple{Array{T,N},Any}} = @abstractmethod

Base.length(::CellArray)::Int = @abstractmethod

function maxsize(self::CellArray{T,N}) where {T,N}
  ms = zeros(Int,N)
  for a in self
    s = size(a)
    @assert length(s) == N
    @inbounds for i in 1:N
      si = ms[i]
      ms[i] = max(si,s[i])
    end
  end
  Tuple(ms)
end

function maxsize(self::CellArray,i::Int)
  s = maxsize(self)
  s[i]
end

maxlength(self::CellArray) = prod(maxsize(self))

Base.eltype(::Type{C}) where C<:CellArray{T,N} where {T,N} = Array{T,N}

function Base.show(io::IO,self::CellArray)
  for (i,a) in enumerate(self)
    println(io,"$i -> $a")
  end
end

function cellsum(self::CellArray{T,N};dim::Int) where {T,N}
  if N == 1
    #TODO This case would lead to a degenerated CellArray{T,0}
    # and deserves special attention
    @notimplemented
  else
    function computesize(a)
      b = [v for (i,v) in enumerate(a) if i!=dim ]
      (b...,)
    end
    function computevals!(a,values)
      #TODO This allocates memory
      sum(a,dims=N)
    end
    CellArrayFromUnaryOp{T,N,T,N-1}(self,computevals!,computesize)
  end
end

function Base.:*(a::CellArray{T,N}, b::CellArray{T,N} ) where {T,N}
  function computevals!(a,b,values)
    if (size(a) != size(values) || size(b) != size(values))
      @notimplemented
    end
    values .= a .* b
    values
  end
  CellArrayFromBinaryOp(a,b,computevals!)
end

"""
Abstract type representing an indexable CellArray.
By implementing a concrete IndexableCellArray, one automatically
gets a type that is also iterable
"""
abstract type IndexableCellArray{T,N} <: CellArray{T,N} end

Base.getindex(::IndexableCellArray{T,N} where {T,N},cell::Int)::Array{T,N} = @abstractmethod

Base.iterate(self::IndexableCellArray) = iterate(self,0)

function Base.iterate(self::IndexableCellArray,state::Int)
  if length(self) == state
    nothing
  else
    k = state+1
    (self[k],k)
  end
end

"""
Abstract type that represents a field with value of type T
evaluated at a collection of points in each cell
"""
const CellFieldValues{T} = CellArray{T,1} where T

# @santiagobadia : Shouldn't we define an object that covers all field types and say that
# T <: FieldType. Idem below and across the project.

"""
Abstract type that represents a function basis with value of type T
evaluated at a collection of points in each cell
"""
const CellBasisValues{T} = CellArray{T,2} where T

"""
An array of values for each cell
"""
const CellValues{T} = CellArray{T,1} where T

"""
An array of points for each cell
"""
const CellPoints{D} = CellValues{Point{D}} where D



# Concrete implementations

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

maxsize(self::ConstantCellArray) = size(self.array)


# TODO: A related abstract type could also be useful
"""
Concrete implementation of `CellArray` that represents the lazy result
of applying a binary operation on two instances of `CellArray`
the functions `computevals!` and `computesize` represents the binary
operation to be performed on the two arrays.
`computesize(::NTuple{P,Int},NTuple{Q,Int})::NTuple{N,Int}`
provides the size of the result and
`computevals!(::Array{A,P},::Array{B,Q},::Array{T,N})::Array{T,N}`
computes the result.
This type is essential for DRY (don't repeat yourself)
"""
struct CellArrayFromBinaryOp{A,P,B,Q,T,N} <: CellArray{T,N}
  a::CellArray{A,P}
  b::CellArray{B,Q}
  computevals!
  computesize
end

function CellArrayFromBinaryOp(a::CellArray{T,N},b::CellArray{T,N},computevals) where {T,N}
  function computesize(a,b)
    @assert a == b
    a
  end
  CellArrayFromBinaryOp{T,N,T,N,T,N}(a,b,computevals,computesize)
end

# @santiagobadia : computesize is returning a CellArray?

function Base.length(self::CellArrayFromBinaryOp)
  @assert length(self.a) == length(self.b)
  length(self.a)
end

function maxsize(self::CellArrayFromBinaryOp{A,P,B,Q,T,N}) where {A,P,B,Q,T,N}
  s = self.computesize(maxsize(self.a), maxsize(self.b))
  @assert isa(s,NTuple{N,Int})
  s
end

function Base.iterate(self::CellArrayFromBinaryOp{A,P,B,Q,T,N}) where {A,P,B,Q,T,N}
  values = Array{T,N}(undef,maxsize(self))
  anext = iterate(self.a)
  bnext = iterate(self.b)
  state = (values,anext,bnext)
  iterate(self,state)
end

function Base.iterate(self::CellArrayFromBinaryOp,state)
  (values,anext,bnext) = state
  if anext == nothing || bnext == nothing
    nothing
  else
    (avals,astate) = anext
    (bvals,bstate) = bnext
    vals = self.computevals!(avals,bvals,values)
    anext = iterate(self.a,astate)
    bnext = iterate(self.b,bstate)
    state = (values,anext,bnext)
    (vals, state)
  end
end

"""
Type that implements the lazy result of an unary operation
on an instance of `CellArray`
"""
struct CellArrayFromUnaryOp{S,M,T,N} <: CellArray{T,N}
  a::CellArray{S,M}
  computevals!
  computesize
end

function CellArrayFromUnaryOp(a::CellArray{T,N},computevals) where {T,N}
  function computesize(a)
    a
  end
  CellArrayFromUnaryOp{T,N,T,N}(a,computevals,computesize)
end

# @santiagobadia : computesize returns a? size(a)?

Base.length(self::CellArrayFromUnaryOp) = length(self.a)

maxsize(self::CellArrayFromUnaryOp) = self.computesize(maxsize(self.a))

function Base.iterate(self::CellArrayFromUnaryOp{S,M,T,N}) where {S,M,T,N}
  values = Array{T,N}(undef,maxsize(self))
  anext = iterate(self.a)
  state = (values,anext)
  iterate(self,state)
end

function Base.iterate(self::CellArrayFromUnaryOp,state)
  (values,anext) = state
  if anext == nothing
    nothing
  else
    (avals,astate) = anext
    vals = self.computevals!(avals,values)
    anext = iterate(self.a,astate)
    state = (values,anext)
    (vals, state)
  end
end
