module CellFields

using Test
using Gridap
using Gridap.Helpers
using Gridap.CellValues: _eq

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
export test_iter_cell_field_without_grad
export test_index_cell_field_without_grad
export test_iter_cell_basis_without_grad
export test_index_cell_basis_without_grad

import Gridap: gradient

const IterCellFieldLike{D,T<:FieldValue,N,R<:FieldLike{D,T,N}} = IterCellValue{R}

const IndexCellFieldLike{D,T<:FieldValue,N,C,R<:FieldLike{D,T,N}} = IndexCellValue{R,C}

const NonIterableCellFieldLike{D,T<:FieldValue,N} = NonIterableCellMap{Point{D},1,T,N}

"""
Umbrella type for CellField and CellBasis
"""
const CellFieldLike{D,T<:FieldValue,N} = Union{IterCellFieldLike{D,T,N},IndexCellFieldLike{D,T,N}}

"""
Returns another CellField or CellBasis object representing the gradient of the given one.
"""
function gradient(::CellFieldLike{D,T,N})::CellFieldLike{D,G,N} where {D,T,G,N}
  @abstractmethod
end

const IterCellField{D,T<:FieldValue,R<:Field{D,T}} = IterCellValue{R}

const IndexCellField{D,T<:FieldValue,C,R<:Field{D,T}} = IndexCellValue{R,C}

"""
Abstract type that represents a cell-wise field, where
`T` stands for the type that represents the field at a point
(e.g., scalar, vector, tensor) and `D` stands for the space
dimension
"""
const CellField{D,T<:FieldValue} = Union{IterCellField{D,T},IndexCellField{D,T}}

function CellField end

const IterCellBasis{D,T<:FieldValue,R<:Basis{D,T}} = IterCellValue{R}

const IndexCellBasis{D,T<:FieldValue,C,R<:Basis{D,T}} = IndexCellValue{R,C}

"""
Abstract type that represents a cell-wise basis for a field space,
where T is the type of value and D the dimension of the domain
"""
const CellBasis{D,T<:FieldValue} = Union{IterCellBasis{D,T},IndexCellBasis{D,T}}

function CellBasis end

"""
Abstract type representing a cellwise transformation between two geometrical
domains
"""
const CellGeomap{D,Z,X,T<:Point{Z,X}} = CellField{D,T}

function CellGeomap end

"""
An array of points for each cell.
This type represents the objects where CellField and CellBasis are evaluated
"""
const CellPoints{D,X} = CellVector{Point{D,X}}

"""
Abstract type that represents a field with value of type T
evaluated at a collection of points in each cell
"""
const CellFieldValues{T<:FieldValue} = CellVector{T}

"""
Abstract type that represents a function basis with value of type T
evaluated at a collection of points in each cell
"""
const CellBasisValues{T<:FieldValue} = CellMatrix{T}

# Testers

function test_iter_cell_field_like(
  f::IterCellFieldLike{D,T,N},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractArray{T,N}},
  g::AbstractArray{<:AbstractArray{G,N}}) where {D,T,G,N}

  _test_iter_cell_field_like(f,x,v,g)

end

function _test_iter_cell_field_like(f,x,v,g)
  test_iter_cell_map(f,x,v)
  fg = gradient(f)
  test_iter_cell_map(fg,x,g)
  _test_field_like_iteration(f,x,v,g)
end

function test_index_cell_field_like(
  f::IndexCellFieldLike{D,T,N},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractArray{T,N}},
  g::AbstractArray{<:AbstractArray{G,N}}) where {D,T,G,N}

  test_index_cell_map(f,x,v)
  fg = gradient(f)
  test_index_cell_map(fg,x,g)
  _test_field_like_iteration(f,x,v,g)
end

function _test_field_like_iteration(m,x,v,g)
  for (mi,xi,vi,gi) in zip(m,x,v,g)
  #  @assert _eq(evaluate(mi,xi),vi)
    @assert _eq(evaluate(gradient(mi),xi),gi)
    @assert typeof(mi) == eltype(m)
  end
end

function test_iter_cell_field(
  f::IterCellField{D,T},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractVector{T}},
  g::AbstractArray{<:AbstractVector{G}}) where {D,T,G}

  _test_iter_cell_field_like(f,x,v,g)
end

function test_iter_cell_field_without_grad(
  f::IterCellField{D,T},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractVector{T}}) where {D,T}

  test_iter_cell_map(f,x,v)
end

function test_index_cell_field(
  f::IndexCellField{D,T},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractVector{T}},
  g::AbstractArray{<:AbstractVector{G}}) where {D,T,G}

  test_index_cell_field_like(f,x,v,g)
end

function test_index_cell_field_without_grad(
  f::IndexCellField{D,T},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractVector{T}}) where {D,T}

  test_index_cell_map(f,x,v)
end

function test_iter_cell_basis(
  f::IterCellBasis{D,T},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractMatrix{T}},
  g::AbstractArray{<:AbstractMatrix{G}}) where {D,T,G}

  _test_iter_cell_field_like(f,x,v,g)
end

function test_iter_cell_basis_without_grad(
  f::IterCellBasis{D,T},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractMatrix{T}}) where {D,T}

  test_iter_cell_map(f,x,v)
end

function test_index_cell_basis(
  f::IndexCellBasis{D,T},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractMatrix{T}},
  g::AbstractArray{<:AbstractMatrix{G}}) where {D,T,G}

  test_index_cell_field_like(f,x,v,g)
end

function test_index_cell_basis_without_grad(
  f::IndexCellBasis{D,T},
  x::CellPoints{D},
  v::AbstractArray{<:AbstractMatrix{T}}) where {D,T}

  test_index_cell_map(f,x,v)
end

end # module
