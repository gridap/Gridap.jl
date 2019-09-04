module DiscFESpaces

using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery

export DiscFESpace

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
Concrete implementation of a discontinuous FE space
with no boundary conditions and for arbitrary
reference elements. Even in the case that the used reference element
has dofs in low dimensional n-faces, the resulting FE space
is always discontinuous.
"""
struct DiscFESpace{D,Z,T} <: FESpace{D,Z,T}
  reffe::RefFE{Z,T}
  nfree::Int
  celldofs::IndexCellVector{Int}
  cellbasis::CellBasis{Z,T}
  trian::Triangulation{D,Z}
end

function DiscFESpace(
  reffe::RefFE{Z,T}, trian::Triangulation{D,Z}) where {D,Z,T}
  cellbasis = _setup_cellbasis(reffe,trian)
  celldofs, nfree = _setup_celldofs(reffe,ncells(trian))
  DiscFESpace{D,Z,T}(reffe,nfree,celldofs,cellbasis,trian)
end

function DiscFESpace(
  reffe::RefFE{D,T}, model::DiscreteModel{D}) where {D,T}

  trian = Triangulation(model)
  DiscFESpace(reffe,trian)
end

num_free_dofs(this::DiscFESpace) = this.nfree

num_diri_dofs(::DiscFESpace) = 0

diri_tags(::DiscFESpace) = Int[]

apply_constraints(
  this::DiscFESpace, cellvec::CellVector, cellids::CellNumber) = cellvec

apply_constraints_rows(
  this::DiscFESpace, cellmat::CellMatrix, cellids::CellNumber) = cellmat

apply_constraints_cols(
  this::DiscFESpace, cellmat::CellMatrix, cellids::CellNumber) = cellmat

function celldofids(this::DiscFESpace)
  this.celldofs
end

function interpolate_values(this::DiscFESpace,fun::Function)
  _setup_interpolate_values(this, fun)
end

function interpolate_diri_values(this::DiscFESpace, funs::Vector{<:Function})
  @assert length(funs) == 0
  T = value_type(this)
  E = dof_type(T)
  zeros(E,0)
end

function CellField(
  this::DiscFESpace{D,Z,T},
  free_dofs::AbstractVector{E},
  diri_dofs::AbstractVector{E}) where {D,Z,T,E}

  @assert length(diri_dofs) == 0

  _setup_cellfield(this,free_dofs)

end

function CellBasis(this::DiscFESpace{D,Z,T}) where {D,Z,T}
  this.cellbasis
end

function Triangulation(this::DiscFESpace{D,Z}) where {Z,D}
  this.trian
end

# Helpers

function _setup_cellbasis(reffe,trian)

  trian_reffes = CellRefFEs(trian)
  @notimplementedif !isa(trian_reffes,ConstantCellValue)
  @assert polytope(reffe).extrusion == polytope(trian_reffes.value).extrusion

  shb = ConstantCellValue(shfbasis(reffe), ncells(trian))
  phi = CellGeomap(trian)
  attachgeomap(shb,phi)
end

function _setup_celldofs(reffe,ncells)
  ncelldofs = length(shfbasis(reffe))
  nfree = ncelldofs*ncells
  celldofs_data = collect(1:nfree)
  celldofs = CellVectorFromDataAndStride(celldofs_data,ncelldofs)
  celldofs, nfree 
end

function _setup_cellfield(this,free_dofs)
  
  ncell = ncells(this.trian)
  ncelldofs = length(shfbasis(this.reffe))
  @assert length(free_dofs) == ncell*ncelldofs

  celldofsvals = CellVectorFromDataAndStride(free_dofs,ncelldofs)

  lincomb(this.cellbasis, celldofsvals)

end

function  _setup_interpolate_values(fesp, fun)

  reffe = fesp.reffe
  dofb = dofbasis(reffe)
  trian = fesp.trian
  phi = CellGeomap(trian)
  uphys = fun ∘ phi

  T = value_type(fesp)
  E = dof_type(T)
  free_dofs = zeros(E, num_free_dofs(fesp))
  diri_dofs = zeros(E, num_diri_dofs(fesp))
  aux = zeros(E, length(dofb))

  _interpolate_kernel!(free_dofs,dofb,uphys,aux)

  (free_dofs, diri_dofs)

end

function  _interpolate_kernel!(free_dofs,dofb,uphys,aux)
  i = 1
  for u in uphys
    evaluate!(dofb,u,aux)
    for a in aux
      free_dofs[i] = a
      i += 1
    end
  end
end

end # module
