module AbstractCellMaps

using Gridap.Helpers
using Gridap.FieldValues
using Gridap.Maps
using Gridap.CellValues

export IterCellMap
export IndexCellMap
export CellMap

export IterCellField
export IndexCellField
export CellField

export IterCellBasis
export IndexCellBasis
export CellBasis

export CellGeomap

export CellFieldValues
export CellBasisValues
export CellPoints

import Gridap: evaluate, gradient
import Gridap: evaluate!, return_size

"""
Abstract object that traverses a set of cells and at every cell returns a
`Map{S,M,T,N}`
"""
const IterCellMap{S,M,T,N,R<:Map{S,M,T,N}} = IterCellValue{R}

"""
Abstract array indexed by cells that returns a `Map{S,M,T,N}`
"""
const IndexCellMap{S,M,T,N,C,R<:Map{S,M,T,N}} = IndexCellValue{R,C}

"""
Abstract object that for a given cell index returns a `Map{S,M,T,N}`
"""
const CellMap{S,M,T,N} = Union{IterCellMap{S,M,T,N},IndexCellMap{S,M,T,N}}

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

"""
Given the maximum size of imput, returns the maximum size of output
"""
function return_size(
  ::CellMap{S,M,T,N},::NTuple{M,Int})::NTuple{N,Int} where {S,M,T,N}
  @abstractmethod
end

"""
Abstract type that represents a cell-wise field, where
`T` stands for the type that represents the field at a point
(e.g., scalar, vector, tensor) and `D` stands for the space
dimension
"""
const IterCellField{D,T,R<:Field{D,T}} = IterCellMap{Point{D},1,T,1,R} where T <:FieldValue

const IndexCellField{D,T,C,R<:Field{D,T}} = IndexCellMap{Point{D},1,T,1,C,R} where T<:FieldValue

const CellField{D,T} = Union{IterCellField{D,T},IndexCellField{D,T}}

function CellField end

"""
Abstract type that represents a cell-wise basis for a field space,
where T is the type of value and D the dimension of the domain
"""
const IterCellBasis{D,T,R<:Basis{D,T}} = IterCellMap{Point{D},1,T,2,R} where T<:FieldValue

const IndexCellBasis{D,T,C,R<:Basis{D,T}} = IndexCellMap{Point{D},1,T,2,C,R} where T<:FieldValue

const CellBasis{D,T} = Union{IterCellBasis{D,T},IndexCellBasis{D,T}}

function CellBasis end

"""
Abstract type representing a cellwise transformation between two geometrical
domains
"""
const CellGeomap{D,Z} = CellField{D,Point{Z}}

function CellGeomap end

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

end # module AbstractCellMaps
