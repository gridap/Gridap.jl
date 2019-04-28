"""
Abstract object that traverses a set of cells and at every cell returns a
`Map{S,M,T,N}`
"""
const IterCellMap{S,M,T,N} = IterCellValue{Map{S,M,T,N}}

"""
Abstract array indexed by cells that returns a `Map{S,M,T,N}`
"""
const IndexCellMap{S,M,T,N,R<:Map{S,M,T,N}} = IndexCellValue{R}

"""
Abstract object that for a given cell index returns a `Map{S,M,T,N}`
"""
const CellMap{S,M,T,N} = Union{IterCellMap{S,M,T,N},IndexCellMap{S,M,T,N}}
# santiagobadia : Problem if IterCellMap and IndexCellMap not same template types?
# Is this correct? IndexCellMap{S,M,T,N} when IndexCellMap{S,M,T,N,R}?

"""
Return the cellwise maps of a `CellMap` on a cellwise set of points
"""
function evaluate(::CellMap{S,M,T,N},::CellArray{S,M})::CellArray{T,N} where {S,M,T,N}
  @abstractmethod
end

"""
Return another `CellMap` object that represents its gradient. Instances of `TG`
have a rank order a unit greater than the ones of `T`
"""
function gradient(::CellMap{S,M,T,N})::CellMap{S,M,TG,N} where {S,M,T<:FieldValue,N,TG}
  @abstractmethod
end

# function Base.show(io::IO,self::CellMap)
#   for (i, a) in enumerate(self)
#     println(io,"$i -> $a")
#   end
# end

"""
Abstract type that represents a cell-wise field, where
`T` stands for the type that represents the field at a point
(e.g., scalar, vector, tensor) and `D` stands for the space
dimension
"""
const IterCellField{D,T} = IterCellMap{Point{D},1,T,1} where {D,T<:FieldValue}
const IndexCellField{D,T,R} = IndexCellMap{Point{D},1,T,1,R} where {D,T<:FieldValue,R}
const CellField{D,T} = Union{IterCellField{D,T},IndexCellField{D,T}}

"""
Abstract type that represents a cell-wise basis for a field space,
where T is the type of value and D the dimension of the domain
"""
const IterCellBasis{D,T} = IterCellMap{Point{D},1,T,2} where {D,T<:FieldValue}
const IndexCellBasis{D,T} = IndexCellMap{Point{D},1,T,2} where {D,T<:FieldValue}
const CellBasis{D,T} = Union{IterCellBasis{D,T},IndexCellBasis{D,T}}

"""
Abstract type representing a cellwise transformation between two geometrical
domains
"""
const CellGeomap{D,Z} = CellField{D,Point{Z}}

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
