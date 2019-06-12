module CellFields

using Test
using Gridap
using Gridap.Helpers

export IterCellFieldLike
export IndexCellFieldLike
export CellFieldLike

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

export test_iter_cell_field
export test_index_cell_field
export test_iter_cell_basis
export test_index_cell_basis

import Gridap: gradient

const IterCellFieldLike = IterCellMap{R} where R<:FieldLike

const IndexCellFieldLike = IndexCellMap{R,D} where {R<:FieldLike,D}

"""
Umbrella type for CellField and CellBasis
"""
const CellFieldLike = Union{IterCellFieldLike{R},IndexCellFieldLike{R,D}} where {R<:FieldLike,D}

"""
Returns another CellField or CellBasis object representing the gradient of the given one.
For efficiency reasons, different cells to this methods should return the same object
"""
function gradient(::CellFieldLike)::CellFieldLike
  @abstractmethod
end

const IterCellField = IterCellMap{R} where R<:Field

const IndexCellField = IndexCellMap{R,D} where {R<:Field,D}

"""
Abstract type that represents a cell-wise field, where
`T` stands for the type that represents the field at a point
(e.g., scalar, vector, tensor) and `D` stands for the space
dimension
"""
const CellField = Union{IterCellField{R},IndexCellField{R,D}} where {R<:Field,D}

function CellField end

const IterCellBasis = IterCellMap{R} where R<:Basis

const IndexCellBasis = IndexCellMap{R,D} where {R<:Basis,D}

"""
Abstract type that represents a cell-wise basis for a field space,
where T is the type of value and D the dimension of the domain
"""
const CellBasis = Union{IterCellBasis{R},IndexCellBasis{R,D}} where {R<:Field,D}

function CellBasis end

"""
Abstract type representing a cellwise transformation between two geometrical
domains
"""
const CellGeomap = CellField{R,C} where {R<:Geomap,C}

function CellGeomap end

"""
An array of points for each cell.
This type represent the objects where CellField and CellBasis are evaluated
"""
const CellPoints{D,X} = CellVector{Point{D,X}} where D

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

# Testers

function test_iter_cell_field_like(
  f::IterCellFieldLike{D,T,N},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractArray{T,N}},
  g::AbstractArray{<:AbstractArray{G,N}}) where {D,T,G,N}

  test_iter_cell_map(f,x,v)
  fg = gradient(f)
  test_iter_cell_map(fg,x,g)
  fg2 = gradient(f)
  @test fg === fg2
end

function test_index_cell_field_like(
  f::IndexCellFieldLike{D,T,N},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractArray{T,N}},
  g::AbstractArray{<:AbstractArray{G,N}}) where {D,T,G,N}

  test_index_cell_map(f,x,v)
  fg = gradient(f)
  test_index_cell_map(fg,x,g)
  fg2 = gradient(f)
  @test fg === fg2
end

function test_iter_cell_field(
  f::IterCellField{D,T},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractVector{T}},
  g::AbstractArray{<:AbstractVector{G}}) where {D,T,G}

  test_iter_cell_field_like(f,x,v,g)
end

function test_index_cell_field(
  f::IndexCellField{D,T},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractVector{T}},
  g::AbstractArray{<:AbstractVector{G}}) where {D,T,G}

  test_index_cell_field_like(f,x,v,g)
end

function test_iter_cell_basis(
  f::IterCellBasis{D,T},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractMatrix{T}},
  g::AbstractArray{<:AbstractMatrix{G}}) where {D,T,G}

  test_iter_cell_field_like(f,x,v,g)
end

function test_index_cell_basis(
  f::IndexCellBasis{D,T},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractMatrix{T}},
  g::AbstractArray{<:AbstractMatrix{G}}) where {D,T,G}

  test_index_cell_field_like(f,x,v,g)
end

end # module
