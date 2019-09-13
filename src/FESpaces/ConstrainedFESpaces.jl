module ConstrainedFESpaces

using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery

export ConstrainedFESpace

import Gridap: num_free_dofs
import Gridap: num_diri_dofs
import Gridap: diri_tags
import Gridap: apply_constraints
import Gridap: apply_constraints_rows
import Gridap: apply_constraints_cols
import Gridap: celldofids
import Gridap: interpolate_values
import Gridap: interpolate_diri_values
import Gridap: CellField
import Gridap: CellBasis
import Gridap: Triangulation

using Base: @propagate_inbounds
import Base: size
import Base: IndexStyle
import Base: getindex

"""
A FESpace that it is constructed from another one plus some constraints
on the free values.
For the moment, we can only fix some free dofs. In the future, we will add
also linear constraints between free dofs.
"""
struct ConstrainedFESpace{D,Z,T,E} <: FESpace{D,Z,T}
  fespace::FESpace{D,Z,T}
  fixeddofs::Vector{Int}
  fixedvals::Vector{E}
  nfree::Int
  celldofs::IndexCellVector{Int}
  is_fixed::Vector{Bool}
  is_free::Vector{Bool}
  dof_to_new_dof::Vector{Int}
end

function ConstrainedFESpace(
  fespace::FESpace{D,Z,T},
  fixeddofs::Vector{Int},
  fixedvals::Vector{E}=zeros(dof_type(T),length(fixeddofs))) where {D,Z,T,E}

  @assert E == dof_type(T)

  celldofs, nfree, is_fixed, is_free, dof_to_new_dof =
    _setup_celldofs(fespace,fixeddofs)
  ConstrainedFESpace{D,Z,T,E}(
    fespace,fixeddofs,fixedvals,nfree,celldofs,
    is_fixed,is_free,dof_to_new_dof)
end

num_free_dofs(this::ConstrainedFESpace) = this.nfree

num_diri_dofs(this::ConstrainedFESpace) = num_diri_dofs(this.fespace)

diri_tags(this::ConstrainedFESpace) = diri_tags(this.fespace)

apply_constraints(
  this::ConstrainedFESpace,
  cellvec::CellVector,
  cellids::CellNumber) = apply_constraints(this.fespace,cellvec,cellids)

apply_constraints_rows(
  this::ConstrainedFESpace,
  cellmat::CellMatrix,
  cellids::CellNumber) = apply_constraints_rows(this.fespace,cellmat,cellids)

apply_constraints_cols(
  this::ConstrainedFESpace,
  cellmat::CellMatrix,
  cellids::CellNumber) = apply_constraints_cols(this.fespace,cellmat,cellids)

function celldofids(this::ConstrainedFESpace)
  this.celldofs
end

function interpolate_values(this::ConstrainedFESpace,fun::Function)
  freevals, dirivals = interpolate_values(this.fespace,fun)
  (freevals[this.is_free], dirivals)
end

function interpolate_values(this::ConstrainedFESpace{D,Z,T},val::T) where {D,Z,T}
  freevals, dirivals = interpolate_values(this.fespace,val)
  (freevals[this.is_free], dirivals)
end

function interpolate_diri_values(this::ConstrainedFESpace, funs::Vector{<:Function})
  interpolate_diri_values(this.fespace,funs)
end

function interpolate_diri_values(this::ConstrainedFESpace{D,Z,T}, vals::Vector{T}) where {D,Z,T}
  interpolate_diri_values(this.fespace,vals)
end

function CellField(
  this::ConstrainedFESpace{D,Z,T},
  free_dofs::AbstractVector{E},
  diri_dofs::AbstractVector{E}) where {D,Z,T,E}

  ndiri = length(diri_dofs)
  v = VectorOfTwoParts(this.dof_to_new_dof,free_dofs,this.fixedvals,ndiri)

  CellField(this.fespace,v,diri_dofs)
end

function CellBasis(this::ConstrainedFESpace{D,Z,T}) where {D,Z,T}
  CellBasis(this.fespace)
end

function Triangulation(this::ConstrainedFESpace{D,Z}) where {Z,D}
  Triangulation(this.fespace)
end

# Helpers

struct VectorOfTwoParts{T,VT1<:AbstractVector{T},VT2<:AbstractVector{T},I} <: AbstractVector{T}
  dof_to_xdof::Vector{I}
  pdof_to_val::VT1
  ndof_to_val::VT2
  offset::I
end

size(v::VectorOfTwoParts) = (length(v.dof_to_xdof),)

IndexStyle(::Type{VectorOfTwoParts{T,I}}) where {T,I} = IndexLinear()

@propagate_inbounds function getindex(v::VectorOfTwoParts,i::Integer)
  @inbounds xdof = v.dof_to_xdof[i]
  if xdof > 0
    @inbounds return v.pdof_to_val[xdof]
  elseif xdof< 0
    @inbounds return v.ndof_to_val[-xdof-v.offset]
  end
end

function _setup_celldofs(fespace,fixeddofs)

  celldofs = celldofids(fespace)
  nfree = num_free_dofs(fespace)
  ndiri = num_diri_dofs(fespace)
  
  nfixed = length(fixeddofs)

  nfree_new = nfree - nfixed

  dof_to_new_dof = zeros(Int,nfree)
  is_fixed = fill(false,nfree)
  is_fixed[fixeddofs] .= true
  is_free = is_fixed.==false

  dof_to_new_dof[is_free] .= 1:nfree_new
  dof_to_new_dof[is_fixed] .= (-(ndiri+1)):-1:(-(ndiri+nfixed))

  celldofs_new = CellVectorFromLocalToGlobalPosAndNeg(
    celldofs,dof_to_new_dof,collect(-1:-1:-ndiri))

  (celldofs_new, nfree_new, is_fixed, is_free, dof_to_new_dof)
end

end # module
