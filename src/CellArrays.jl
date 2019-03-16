export CellArray, IndexableCellArray
export CellFieldValues, CellBasisValues
export CellScalars, CellVectors, CellMatrices

"""
Abstract type representing an iterable collection of Arrays{T,N},
where each array is associated with a cell.
"""
abstract type CellArray{T,N} end

Base.iterate(::CellArray)::Union{Nothing,Tuple{Array{T,N},Any}} = @abstractmethod

Base.iterate(::CellArray,state)::Union{Nothing,Tuple{Array{T,N},Any}} = @abstractmethod

Base.length(::CellArray)::Int = @abstractmethod

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

Base.length(::IndexableCellArray)::Int = @abstractmethod

Base.iterate(self::IndexableCellArray) = iterate(self,0)

function Base.iterate(self::IndexableCellArray,state::Int)
  if length(self) == state
    nothing
  else
    k = state+1
    (self[k],k)
  end
end

# TODO: To discuss. I (@fverdugo) think that it is not required to
# define new abstract types to represent the following concepts. They can be
# expressed all of them as CellArrays for specific values of the N parameter

"""
Abstract type that represents a field with value of type T
evaluated at a collection of points in each cell
"""
const CellFieldValues{T} = CellArray{T,1} where T

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

