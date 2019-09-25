module ZeroMeanFESpaces

using Gridap
using Gridap.Helpers
using Gridap.ConstrainedFESpaces: VectorOfTwoParts

export ZeroMeanFESpace

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

struct ZeroMeanFESpace{D,Z,T} <: FESpace{D,Z,T}
  fespace::ConstrainedFESpace{D,Z,T}
  voli::Vector{Float64}
  vol::Float64
  quadorder::Int
end

function ZeroMeanFESpace(fespace::FESpace{D,Z,T},quadorder::Int) where {D,Z,T}

  fixeddofs = [1,]
  _fespace = ConstrainedFESpace(fespace,fixeddofs)

  vol_i, vol = _setup_vols(fespace,quadorder)

  ZeroMeanFESpace(_fespace,vol_i,vol,quadorder)

end

num_free_dofs(this::ZeroMeanFESpace) = num_free_dofs(this.fespace)

num_diri_dofs(this::ZeroMeanFESpace) = num_diri_dofs(this.fespace)

diri_tags(this::ZeroMeanFESpace) = diri_tags(this.fespace)

apply_constraints(
  this::ZeroMeanFESpace,
  cellvec::CellVector,
  cellids::CellNumber) = apply_constraints(this.fespace,cellvec,cellids)

apply_constraints_rows(
  this::ZeroMeanFESpace,
  cellmat::CellMatrix,
  cellids::CellNumber) = apply_constraints_rows(this.fespace,cellmat,cellids)

apply_constraints_cols(
  this::ZeroMeanFESpace,
  cellmat::CellMatrix,
  cellids::CellNumber) = apply_constraints_cols(this.fespace,cellmat,cellids)

function celldofids(this::ZeroMeanFESpace)
  celldofids(this.fespace)
end

function interpolate_values(this::ZeroMeanFESpace,fun::Function)
  interpolate_values(this.fespace,fun)
end

function interpolate_values(this::ZeroMeanFESpace{D,Z,T},val::T) where {D,Z,T}
    interpolate_values(this.fespace,val)
end

function interpolate_diri_values(this::ZeroMeanFESpace, funs::Vector{<:Function})
  interpolate_diri_values(this.fespace,funs)
end

function interpolate_diri_values(this::ZeroMeanFESpace{D,Z,T}, vals::Vector{T}) where {D,Z,T}
  interpolate_diri_values(this.fespace,vals)
end

function CellField(
  this::ZeroMeanFESpace{D,Z,T},
  free_dofs::AbstractVector{E},
  diri_dofs::AbstractVector{E}) where {D,Z,T,E}

  ndiri = length(diri_dofs)
  v = VectorOfTwoParts(
    this.fespace.dof_to_new_dof,free_dofs,this.fespace.fixedvals,ndiri)

  fixedval = _compute_new_fixedval(v,this.voli,this.vol)
  fixedvals = [fixedval,]

  free_dofs[:] .+= fixedval

  w = VectorOfTwoParts(
    this.fespace.dof_to_new_dof,free_dofs,fixedvals,ndiri)

  CellField(this.fespace.fespace,w,diri_dofs)
end

function CellBasis(this::ZeroMeanFESpace{D,Z,T}) where {D,Z,T}
  CellBasis(this.fespace)
end

function Triangulation(this::ZeroMeanFESpace{D,Z}) where {Z,D}
  Triangulation(this.fespace)
end

# Helpers

function _setup_vols(fespace,quadorder)

  @notimplementedif length(diri_tags(fespace)) != 0
  @notimplementedif ! (value_type(fespace) <: Real)

  V = TestFESpace(fespace)
  U = TrialFESpace(fespace)

  assem = SparseMatrixAssembler(V,U)

  trian = Triangulation(fespace)
  quad = CellQuadrature(trian,degree=quadorder)

  bh = FEBasis(V)

  cv = integrate(bh.cellbasis,trian,quad)
  cn = IdentityCellNumber(Int,ncells(trian))

  vol_i = assemble(assem,(cv,cn))

  vol = sum(vol_i)

  vol_i, vol 
end

function _compute_new_fixedval(v,vol_i,vol)

  c = 0.0

  for (i,vi) in enumerate(v)
    c += vi*vol_i[i]
  end

  c = -c/vol

  c

end

end # module
