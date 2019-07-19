module FESpaces

using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery

export free_dofs
export diri_dofs
export FESpace
export TestFESpace
export TrialFESpace
export num_free_dofs
export num_diri_dofs
export diri_tags
export apply_constraints
export apply_constraints_rows
export apply_constraints_cols
export celldofids
export interpolated_values
export interpolated_diri_values
export value_type

export ConformingFESpace
export FESpaceWithDirichletData

import Gridap: CellField
import Gridap: CellBasis

"""
Abstract FE Space parameterized with respec to the environment dimension `D`,
the cell dimension `Z`, and the field type `T` (rank and entry type).
A FE space has to be understood as the direct sum of two spaces, the one with
arbitrary free values and zero Dirichlet data and the one with zero free values
and a given Dirichlet data
"""
abstract type FESpace{D,Z,T} end

num_free_dofs(::FESpace)::Int = @abstractmethod

num_diri_dofs(::FESpace)::Int = @abstractmethod

diri_tags(::FESpace)::Vector{Int} = @abstractmethod

function apply_constraints(::FESpace, cellvec::CellVector, cellids::CellNumber)::CellVector
  @abstractmethod
end

function apply_constraints_rows(::FESpace, cellmat::CellMatrix, cellids::CellNumber)::CellMatrix
  @abstractmethod
end

function apply_constraints_cols(::FESpace, cellmat::CellMatrix, cellids::CellNumber)::CellMatrix
  @abstractmethod
end

"""
Cell DOFs ids after applying constraints
"""
function celldofids(::FESpace)::CellVector{Int}
  @abstractmethod
end

"""
Interpolates a function and returns a vector of free values and another vector
of dirichlet ones. It is templatized with respect to E, the type of field
value in the `FESpace{D,Z,T}` at hand, i.e.,  `E = eltype(T)`
"""
function interpolated_values(::FESpace,::Function)::Tuple{Vector{E},Vector{E}} where E
  @abstractmethod
end
# @santiagobadia : Private method
# @santiagobadia : interpolated -> interpolate

"""
Given a FESpace and its array of labels that represent the Dirichlet boundary
in the geometrical model and an array of functions for every Dirichlet
boundary label, this method interpolates all these functions on their
boundaries, and provides the global vector of Dirichlet DOFs
"""
function interpolated_diri_values(::FESpace, funs::Vector{<:Function})::Vector{E} where E
  @abstractmethod
end
# @santiagobadia : Private method
# @santiagobadia : interpolated -> interpolate

function interpolated_diri_values(this::FESpace, fun::Function)
  tags = diri_tags(this)
  interpolated_diri_values(this,fill(fun,length(tags)))
end
# @santiagobadia : Private method
# @santiagobadia : interpolated -> interpolate

"""
Returns the CellField that represents the FE function thorough its free and
dirichlet values. E = eltype(T)
"""
function CellField(
  ::FESpace{D,Z,T},free_dofs::AbstractVector{E},diri_dofs::AbstractVector{E})::CellField{Z,T} where {D,Z,T,E}
  @abstractmethod
end
# @santiagobadia : Private method (?), the user will use FEFunction constructor

function CellBasis(::FESpace{D,Z,T})::CellBasis{Z,T} where {D,Z,T}
  @abstractmethod
end
# @santiagobadia : Private method

value_type(::FESpace{D,Z,T}) where {D,Z,T} = T
# @santiagobadia : Private method

function TestFESpace(this::FESpace{D,Z,T}) where {D,Z,T}
  E = eltype(T)
  dv = zeros(E,num_diri_dofs(this))
  return FESpaceWithDirichletData(this, dv)
end

function TrialFESpace( this::FESpace, funs::Vector{<:Function}) where {D}
  tags = diri_tags(this)
  @assert length(tags) == length(funs)
  dv = interpolated_diri_values(this,funs)
  return FESpaceWithDirichletData(this, dv)
end

function TrialFESpace( this::FESpace, fun::Function) where {D}
  dv = interpolated_diri_values(this,fun)
  return FESpaceWithDirichletData(this, dv)
end

"""
FESpace whose Dirichlet component has been constrained
"""
struct FESpaceWithDirichletData{D,Z,T,E,V<:FESpace{D,Z,T}} <: FESpace{D,Z,T}
  fespace::V
  diri_dofs::Vector{E}
end

function free_dofs end

diri_dofs(f::FESpaceWithDirichletData) = f.diri_dofs

num_free_dofs(f::FESpaceWithDirichletData) = num_free_dofs(f.fespace)

num_diri_dofs(f::FESpaceWithDirichletData) = num_diri_dofs(f.fespace)

diri_tags(f::FESpaceWithDirichletData) = diri_tags(f.fespace)

function apply_constraints(
  f::FESpaceWithDirichletData, cellvec::CellVector, cellids::CellNumber)
  apply_constraints(f.fespace,cellvec,cellids)
end

function apply_constraints_rows(
  f::FESpaceWithDirichletData, cellmat::CellMatrix, cellids::CellNumber)
  apply_constraints_rows(f.fespace,cellmat,cellids)
end

function apply_constraints_cols(
  f::FESpaceWithDirichletData, cellmat::CellMatrix, cellids::CellNumber)
  apply_constraints_cols(f.fespace,cellmat,cellids)
end

function celldofids(f::FESpaceWithDirichletData)
  celldofids(f.fespace)
end

function interpolated_values(f::FESpaceWithDirichletData,fun::Function)
  free_vals, _ = interpolated_values(f.fespace,fun)
  free_vals, f.diri_dofs
end

function CellField(
  f::FESpaceWithDirichletData{D,Z,T},free_dofs::AbstractVector{E},diri_dofs::AbstractVector{E})where {D,Z,T,E}
  CellField(f.fespace,free_dofs,f.diri_dofs)
end

function CellBasis(f::FESpaceWithDirichletData)
  CellBasis(f.fespace)
end

function interpolated_diri_values(this::FESpaceWithDirichletData, funs::Vector{<:Function})
  f.diri_dofs
end

"""
Conforming FE Space, where only one RefFE is possible in the whole mesh
"""
struct ConformingFESpace{D,Z,T} <: FESpace{D,Z,T}
  dim_to_nface_eqclass::Vector{<:IndexCellArray{Int}}
  cell_eqclass::IndexCellArray{Int}
  num_free_dofs::Int
  num_diri_dofs::Int
  diri_tags::Vector{Int}
  _reffes::LagrangianRefFE{D,T}
  _triangulation::Triangulation{D,Z}
  _gridgraph::FullGridGraph
  _labels::FaceLabels
  _basis::CellBasis{Z,T}
end

function ConformingFESpace(
  reffe::LagrangianRefFE{D,T},
  trian::Triangulation{D,Z},
  graph::FullGridGraph,
  labels::FaceLabels,
  diri_tags::Vector{Int}) where {D,Z,T}
  args = _setup_conforming_fe_fields(reffe,trian,graph,labels,diri_tags,D)
  ConformingFESpace{D,Z,T}(args...)
end

function ConformingFESpace(
  reffe::LagrangianRefFE{D,T},
  trian::Triangulation{D,Z},
  graph::FullGridGraph,
  labels::FaceLabels) where {D,Z,T}
  return ConformingFESpace(reffe, trian, graph, labels, ())
end

function ConformingFESpace(::Type{T},model::DiscreteModel{D},order,diri_tags) where {D,T}
  grid = Grid(model,D)
  trian = Triangulation(grid)
  graph = FullGridGraph(model)
  labels = FaceLabels(model)
  orders = fill(order,D)
  polytope = _polytope(celltypes(grid))
  fe = LagrangianRefFE{D,T}(polytope, orders)
  _diri_tags = _setup_diri_tags(model,diri_tags)
  ConformingFESpace(fe,trian,graph,labels,_diri_tags)
end

num_free_dofs(this::ConformingFESpace) = this.num_free_dofs

num_diri_dofs(this::ConformingFESpace) = this.num_diri_dofs

diri_tags(f::ConformingFESpace) = f.diri_tags

function apply_constraints(
  this::ConformingFESpace, cellvec::CellVector, cellids::CellNumber)
  cellvec
end

function apply_constraints_rows(
  this::ConformingFESpace, cellmat::CellMatrix, cellids::CellNumber)
  cellmat
end

function apply_constraints_cols(
  this::ConformingFESpace, cellmat::CellMatrix, cellids::CellNumber)
  cellmat
end

function celldofids(this::ConformingFESpace)
  this.cell_eqclass
end

function interpolated_values(this::ConformingFESpace,f::Function)
  _interpolated_values(this,f)
end

function interpolated_diri_values(this::ConformingFESpace, funs::Vector{<:Function})
  _interpolated_diri_values(this,funs)
end

function CellField(
  fespace::ConformingFESpace{D,Z,T},
  free_dofs::AbstractVector{E},
  diri_dofs::AbstractVector{E}) where {D,Z,T,E}

  @assert E == eltype(T)
  @assert num_free_dofs(fespace) == length(free_dofs)
  @assert num_diri_dofs(fespace) == length(diri_dofs)
  reffe = fespace._reffes
  trian = fespace._triangulation
  celldofs = fespace.cell_eqclass
  shb = CellBasis(fespace)
  cdofs = CellVectorFromLocalToGlobalPosAndNeg(celldofs, free_dofs, diri_dofs)
  lincomb(shb,cdofs)
end

CellBasis(this::ConformingFESpace) = this._basis

# Helpers

_setup_diri_tags(model,tags) = tags

function _setup_diri_tags(model,name::String)
  _setup_diri_tags(model,[name,])
end

function _setup_diri_tags(model,names::Vector{String})
  [ tag_from_name(model,s) for s in names ]
end

_polytope(celltypes) = @notimplemented

function _polytope(celltypes::ConstantCellValue)
  code = celltypes.value
  Polytope(code)
end

function _interpolated_values(fesp::ConformingFESpace{D,Z,T},fun::Function) where {D,Z,T}
  reffe = fesp._reffes
  dofb = reffe.dofbasis
  trian = fesp._triangulation
  phi = CellGeomap(trian)
  uphys = fun ∘ phi
  celldofs = fesp.cell_eqclass
  nfdofs = fesp.dim_to_nface_eqclass
  maxs = max([length(nfdofs[i]) for i=1:D+1]...)
  E = eltype(T)
  free_dofs = zeros(E, num_free_dofs(fesp))
  diri_dofs = zeros(E, num_diri_dofs(fesp))
  aux = zeros(E, maxs)
  for (imap,l2g) in zip(uphys,celldofs)
    evaluate!(dofb,imap,aux)
    for (i,gdof) in enumerate(l2g)
      if (gdof > 0)
        free_dofs[gdof] = aux[i]
      else
        diri_dofs[-gdof] = aux[i]
      end
    end
  end
  return free_dofs, diri_dofs
end

function _interpolated_diri_values(fesp::ConformingFESpace{D,Z,T},funs) where {D,Z,T}
  labels = fesp._labels
  nf_labs_all = [ labels_on_dim(labels,idim) for idim in 0:D]
  nf_dofs_all = fesp.dim_to_nface_eqclass
  dtags = fesp.diri_tags
  @assert length(dtags) == length(funs)
  E = eltype(T)
  diri_dofs_all = zeros(E, num_diri_dofs(fesp))
  for (ifunc,f) in enumerate(funs)
    _ , fh_diri_dofs = interpolated_values(fesp,f)
    # @santiagobadia: For performance issues, implement a new interpolate
    # restricted to cells on the boundary for performance
    for idim in 0:D
      nf_labs = nf_labs_all[idim+1]
      nf_dofs = nf_dofs_all[idim+1]
      # How to extract this part? Do it correctly, with a FEFunction
      for (nf,nflab) in enumerate(nf_labs)
        if (_is_diri(nflab, (dtags[ifunc],), labels))
          for dof in nf_dofs[nf]
            dof *= -1
            diri_dofs_all[dof] = fh_diri_dofs[dof]
          end
        end
      end
    end
  end
  return diri_dofs_all
end

function _setup_conforming_fe_fields(reffe,trian,graph,labels,diri_tags,D)
  dim_to_nface_eqclass, nfree, ndiri  = _generate_dim_to_nface_eqclass(
    reffe, graph, labels, diri_tags)
  cellvefs_dim = [connections(graph,D,i) for i in 0:D]
  offset = length.(dim_to_nface_eqclass)
  for i in 2:length(offset)
    offset[i] += offset[i-1]
  end
  offset = tuple(offset[1:(end-1)]...)
  cellvefs = local_append(offset, cellvefs_dim...)
  dofs_all = append(dim_to_nface_eqclass...)
  cell_eqclass = CellVectorByComposition(cellvefs, dofs_all)
  shb = ConstantCellValue(reffe.shfbasis, ncells(trian))
  phi = CellGeomap(trian)
  basis = attachgeomap(shb,phi)
  return dim_to_nface_eqclass, cell_eqclass, nfree, ndiri, diri_tags,
    reffe, trian, graph, labels, basis
end

function _generate_dim_to_nface_eqclass(
  reffe::RefFE{D,T},
  graph::FullGridGraph,
  labels::FaceLabels,
  diri_tags::Vector{Int}) where {D,T}

  dim_to_nface_eqclass = Int[]
  c=1
  c_n = -1
  for vef_dim in 0:D
    vefcells= connections(graph,vef_dim,D)
    cellvefs= connections(graph,D,vef_dim)
    vef_labels = labels_on_dim(labels,vef_dim)
    num_vefs = length(vefcells)
    nfdofs=Array{Array{Int64},1}(undef,num_vefs)
    nfdofs_l = Int[]
    nfdofs_g = zeros(Int, num_vefs+1)
    nfdofs_g[1] = 1
    for (ignf,nf) in enumerate(vefcells)
      owner_cell = nf[1]
      lid_vef_dim = findfirst(i->i==ignf,cellvefs[owner_cell])
      lid_vef = reffe.polytope.nf_dim[end][vef_dim+1][1]+lid_vef_dim-1
      # @santiagobadia : Better a method for nfs of a particular type...
      num_nf_dofs = length(reffe.nfacedofs[lid_vef])
      if ( _is_diri(vef_labels[ignf],diri_tags,labels) )
        nfdofs_l = Int[nfdofs_l..., c_n:-1:c_n-num_nf_dofs+1... ]
        c_n -= num_nf_dofs
      else
        nfdofs_l = Int[nfdofs_l..., c:c+num_nf_dofs-1... ]
        c += num_nf_dofs
      end
      nfdofs_g[ignf+1] += num_nf_dofs + nfdofs_g[ignf]
    end
    dim_to_nface_eqclass = [
      dim_to_nface_eqclass..., CellVectorFromDataAndPtrs(nfdofs_l, nfdofs_g) ]
  end
  return [ dim_to_nface_eqclass , c-1, -c_n-1 ]
end

@inline function _is_diri(v,dt,labels)
  for tag in dt
    labs = labels_on_tag(labels,tag)
    for label in labs
      if (label == v)
        return true
      end
    end
  end
  return false
end

end # module FESpaces
