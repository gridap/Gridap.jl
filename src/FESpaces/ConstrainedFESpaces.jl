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
end

function ConstrainedFESpace(
  fespace::FESpace{D,Z,T},
  fixeddofs::Vector{Int},
  fixedvals::Vector{E}=zeros(dof_type(T),length(fixeddofs))) where {D,Z,T,E}

  @assert E == dof_type(T)

  celldofs, nfree, is_fixed, is_free = _setup_celldofs(fespace,fixeddofs)
  ConstrainedFESpace{D,Z,T,E}(
    fespace,fixeddofs,fixedvals,nfree,celldofs,is_fixed,is_free)
end

num_free_dofs(this::ConstrainedFESpace) = this.nfree

num_diri_dofs(this::ConstrainedFESpace) = num_diri_dofs(this.fespace)

diri_tags(this::ConstrainedFESpace) = diri_tags(this.fespace)

apply_constraints(
  this::ConstrainedFESpace,
  cellvec::CellVector,
  cellids::CellNumber) = cellvec

apply_constraints_rows(
  this::ConstrainedFESpace,
  cellmat::CellMatrix,
  cellids::CellNumber) = cellmat

apply_constraints_cols(
  this::ConstrainedFESpace,
  cellmat::CellMatrix,
  cellids::CellNumber) = cellmat

function celldofids(this::ConstrainedFESpace)
  this.celldofs
end

function interpolate_values(this::ConstrainedFESpace,fun::Function)
  freevals, dirivals = interpolate_values(this.fespace,fun)
  (freevals[this.is_free], dirivals)
end

function interpolate_diri_values(this::ConstrainedFESpace, funs::Vector{<:Function})
  interpolate_diri_values(this.fespace,funs)
end

function CellField(
  this::ConstrainedFESpace{D,Z,T},
  free_dofs::AbstractVector{E},
  diri_dofs::AbstractVector{E}) where {D,Z,T,E}

  nfree_old = num_free_dofs(this.fespace)

  free_dofs_old = zeros(E,nfree_old)
  free_dofs_old[this.is_free] .= free_dofs
  free_dofs_old[this.is_fixed] .= this.fixedvals

  CellField(this.fespace,free_dofs_old,diri_dofs)
end

function CellBasis(this::ConstrainedFESpace{D,Z,T}) where {D,Z,T}
  CellBasis(this.fespace)
end

function Triangulation(this::ConstrainedFESpace{D,Z}) where {Z,D}
  Triangulation(this.fespace)
end

# Helpers

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

  (celldofs_new, nfree_new, is_fixed, is_free)
end

end # module
