
"""
An abstract type that maps a CellArray into another CellArray
"""
abstract type CellFunction{S,M,T,N} end

"""
Returns the evaluation of a `CellFunction`
"""
function evaluate(::CellFunction{S,M,T,N},::CellArray{S,M})::CellArray{T,N} where {S,M,T,N}
  @abstractmethod
end

"""
Returns another CellFunction object that represents its gradient.
Instances of `TG` have a rank order a unit greater than the ones of `T`
"""
function gradient(::CellFunction{S,M,T,N})::CellFunction{S,M,TG,N} where {S,M,T<:FieldValue,N,TG}
  @abstractmethod
end

"""
Abstract type that represents a cell-wise basis for a field space,
where T is the type of value and D the dimension of the domain
"""
const CellBasis{D,T} = CellFunction{Point{D},1,T,2} where {D,T<:FieldValue}

"""
Abstract type that represents a cell-wise field, where
`T` stands for the type that represents the field at a point
(e.g., scalar, vector, tensor) and `D` stands for the space
dimension
"""
const CellField{D,T} = CellFunction{Point{D},1,T,1} where {D,T<:FieldValue}

"""
Abstract type representing a cell-wise transformation
between two geometrical domains
"""
const CellGeomap{D,Z} = CellField{D,Point{Z}}

# Abstract types for the input and output values
# of CellFields and CellBasis

"""
An array of points for each cell.
This type represent the objects where CellField and CellBasis are evaluated
"""
const CellPoints{D} = CellVector{Point{D}} where D

"""
Abstract type that represents a field with value of type T
evaluated at a collection of points in each cell
"""
const CellFieldValues{T} = CellVector{T} where T <: FieldValue

"""
Abstract type that represents a function basis with value of type T
evaluated at a collection of points in each cell
"""
const CellBasisValues{T} = CellArray{T,2} where T <: FieldValue
