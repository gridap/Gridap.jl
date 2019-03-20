export CellArray, IndexableCellArray
export maxlength, maxsize
export CellFieldValues, CellBasisValues
export CellPoints
export CellScalars, CellVectors, CellMatrices
export ConstantCellArray

"""
Abstract type representing an iterable collection of Arrays{T,N},
where each array is associated with a cell.
"""
abstract type CellArray{T,N} end

Base.iterate(::CellArray)::Union{Nothing,Tuple{Array{T,N},Any}} = @abstractmethod

Base.iterate(::CellArray,state)::Union{Nothing,Tuple{Array{T,N},Any}} = @abstractmethod

Base.length(::CellArray)::Int = @abstractmethod

function maxsize(self::CellArray{N}) where N
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

maxlength(self::CellArray) = prod(maxsize(self))

Base.eltype(::Type{C}) where C<:CellArray{T,N} where {T,N} = Array{T,N}

function Base.show(io::IO,self::CellArray)
  for (i,a) in enumerate(self)
    println(io,"$i -> $a")
  end
end

"""
Abstract type representing an indexable CellArray.
By implementing a concrete IndexableCellArray, one automatically
gets an type that is also iterable
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

"""
An array of points for each cell
"""
const CellPoints{D} = CellFieldValues{Point{D}} where D

"""
Abstract type that represents a function basis with value of type T
evaluated at a collection of points in each cell
"""
const CellBasisValues{T} = CellArray{T,2} where T

"""
Abstract type that represents a scalar of value T
associated with a collection of points in each cell
"""
const CellScalars{T} = CellArray{T,1} where T

"""
Abstract type that represents a vector of value T
associated with a collection of points in each cell
(typically the cell rhs vector at the quadrature points)
"""
const CellVectors{T} = CellArray{T,2} where T

"""
Abstract type that represents a matrix of value T
associated with a collection of points in each cell
(typically the cell matrix at the quadrature points)
"""
const CellMatrices{T} = CellArray{T,3} where T


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
