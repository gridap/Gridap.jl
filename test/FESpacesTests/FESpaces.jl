module FESpaces

using Gridap
using Gridap.Helpers
using Gridap.RefFEs
using Gridap.Polytopes
using Gridap.CellValues
using Gridap.CellMaps
using Gridap.CellValues.ConstantCellValues
using Gridap.CellValues.Wrappers
using Gridap.CellValues.Append
using Gridap.Geometry
using Gridap.FieldValues
using Gridap.CellMaps.Operations: CellFieldFromExpand

export FEFunction
export free_dofs
export diri_dofs
export FESpace
export num_free_dofs
export num_diri_dofs
export apply_constraints
export apply_constraints_rows
export apply_constraints_cols
export interpolated_values

import Gridap: evaluate, gradient, return_size

export ConformingFESpace
export ConformingFEFunction

"""
Abstract type representing a FE Function
A FE function is a member of a FESpace.
Since a FE space is the direct sum of free + dirichlet spaces
a FE function is uniquely represented by two vectors of components aka
`free_dofs` and `diri_dofs` which are Vectors{E} with
E = eltype(T)
"""
abstract type FEFunction{Z,T,R} <: IterCellField{Z,T,R} end

function free_dofs(::FEFunction)::Vector{E} where E
  @abstractmethod
end

function diri_dofs(::FEFunction)::Vector{E} where E
  @abstractmethod
end

"""
Abstract FE Space parameterized with respec to the environment dimension `D`,
the cell dimension `Z`, and the field type `T` (rank)
A FE space has to be understood as the direct sum of two spaces free + dirichlet
"""
abstract type FESpace{D,Z,T} end

num_free_dofs(::FESpace)::Int = @abstractmethod

num_diri_dofs(::FESpace)::Int = @abstractmethod

function apply_constraints(
  ::FESpace, cellvec::CellVector)::Tuple{CellVector,CellVector{Int}}
  @abstractmethod
end

function apply_constraints_rows(
  ::FESpace, cellmat::CellMatrix)::Tuple{CellMatrix,CellVector{Int}}
  @abstractmethod
end

function apply_constraints_cols(
  ::FESpace, cellmat::CellMatrix)::Tuple{CellMatrix,CellVector{Int}}
  @abstractmethod
end

"""
Interpolates a function and returns a vector of free values and another of dirichlet ones
E = eltype(T)
"""
function interpolated_values(::FESpace,::Function)::Tuple{Vector{E},Vector{E}} where E
  @abstractmethod
end

"""
Returns the FE function represented be the  free and dirichlet values
E = eltype(T)
"""
function FEFunction(::FESpace,free_dofs::Vector{E},diri_dofs::Vector{E})::FEFunction where E
  @abstractmethod
end

function interpolate(this::FESpace,fun::Function)
  free_vals, diri_vals = interpolated_values(this,fun)
  FEFunction(this,free_vals,diri_vals)
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
end

function ConformingFESpace(
  reffe::LagrangianRefFE{D,T},
  trian::Triangulation{D,Z},
  graph::FullGridGraph,
  labels::FaceLabels) where {D,Z,T}
  return ConformingFESpace(reffe, trian, graph, labels, ())
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

num_free_dofs(this::ConformingFESpace) = this.num_free_dofs

num_diri_dofs(this::ConformingFESpace) = this.num_diri_dofs

function apply_constraints(
  this::ConformingFESpace, cellvec::CellVector)
  return cellvec, this.cell_eqclass
end

function apply_constraints_rows(
  this::ConformingFESpace, cellmat::CellMatrix)
  return cellmat, this.cell_eqclass
end

function apply_constraints_cols(
  this::ConformingFESpace, cellmat::CellMatrix)
  return cellmat, this.cell_eqclass
end

function interpolated_values(this::ConformingFESpace,f::Function)
  _interpolated_values(this,f)
end

function FEFunction(
  this::ConformingFESpace,free_dofs::Vector{E},diri_dofs::Vector{E}) where E
  ConformingFEFunction(this,free_dofs,diri_dofs)
end

function interpolate(this::ConformingFESpace,f::Function)
  _interpolate(this,f)
end

struct ConformingFEFunction{
   Z,T,R,
   F<:CellFieldFromExpand{
    Z,T,<:CellBasis{Z},<:CellVectorFromLocalToGlobalPosAndNeg,R
   }} <: FEFunction{Z,T,R}
  cfield::F
end

function ConformingFEFunction(
  fespace::ConformingFESpace{D,Z,T},
  free_dofs::Vector{E},
  diri_dofs::Vector{E}) where {D,Z,T,E}

  @assert E == eltype(T)
  @assert num_free_dofs(fespace) == length(free_dofs)
  @assert num_diri_dofs(fespace) == length(diri_dofs)
  reffe = fespace._reffes
  trian = fespace._triangulation
  celldofs = fespace.cell_eqclass
  shb = ConstantCellValue(reffe.shfbasis, ncells(trian))
  cdofs = CellVectorFromLocalToGlobalPosAndNeg(celldofs, free_dofs, diri_dofs)
  cfield = CellFieldFromExpand(shb, cdofs)
  ConformingFEFunction(cfield)
end

function free_dofs(this::ConformingFEFunction)
  this.cfield.coeffs.gid_to_val_pos.v
end

function diri_dofs(this::ConformingFEFunction)
  this.cfield.coeffs.gid_to_val_neg.v
end

function evaluate(this::ConformingFEFunction{Z},q::CellPoints{Z}) where Z
  evaluate(this.cfield,q)
end

# Helpers

function _interpolated_values(fesp::ConformingFESpace{D,Z,T},fun::Function) where {D,Z,T}
  reffe = fesp._reffes
  dofb = reffe.dofbasis
  trian = fesp._triangulation
  phi = geomap(trian)
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

function _setup_conforming_fe_fields(reffe,trian,graph,labels,diri_tags,D)
  dim_to_nface_eqclass, nfree, ndiri  = _generate_dim_to_nface_eqclass(
    reffe, graph, labels, diri_tags)
  cellvefs_dim = [connections(graph,D,i) for i in 0:D]
  offset = length.(dim_to_nface_eqclass)
  for i in 2:length(offset)
    offset[i] += offset[i-1]
  end
  offset = tuple(offset...)
  cellvefs = IndexCellValueByLocalAppendWithOffset(offset, cellvefs_dim...)
  dofs_all = IndexCellValueByGlobalAppend(dim_to_nface_eqclass...)
  cell_eqclass = CellVectorByComposition(cellvefs, dofs_all)
  return dim_to_nface_eqclass, cell_eqclass, nfree, ndiri, diri_tags,
    reffe, trian, graph, labels
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
