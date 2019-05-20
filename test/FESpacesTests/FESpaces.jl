module FESpaces

using Gridap
using Gridap.Helpers
using Gridap.RefFEs
using Gridap.Polytopes
using Gridap.CellValues
using Gridap.CellValues.ConstantCellValues
using Gridap.CellValues.Wrappers
using Gridap.CellValues.Append
using Gridap.Geometry
using Gridap.FieldValues

export FESpace
export num_free_dofs
export num_fixed_dofs
export apply_constraints
export apply_constraints_rows
export apply_constraints_cols

export ConformingFESpace

"""
Abstract FE Space parameterized with respec to the environment dimension `D`,
the cell dimension `Z`, and the field type `T` (rank)
"""
abstract type FESpace{D,Z,T} end

num_free_dofs(::FESpace)::Int = @abstractmethod

num_fixed_dofs(::FESpace)::Int = @abstractmethod

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
Conforming FE Space, where only one RefFE is possible in the whole mesh
"""
struct ConformingFESpace{D,Z,T} <: FESpace{D,Z,T}
  dim_to_nface_eqclass::Vector{<:IndexCellArray{Int}}
  cell_eqclass::IndexCellArray{Int}
  num_free_dofs::Int
  num_fixed_dofs::Int
  dir_tags::Vector{Int}
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
  dir_tags::Vector{Int}) where {D,Z,T}

  dim_to_nface_eqclass, nfree, nfixed  = _generate_dim_to_nface_eqclass(
    reffe, graph, labels, dir_tags)

  cellvefs_dim = [connections(graph,D,i) for i in 0:D]
  offset = length.(dim_to_nface_eqclass)
  for i in 2:length(offset)
    offset[i] += offset[i-1]
  end
  offset = tuple(offset...)
  cellvefs = IndexCellValueByLocalAppendWithOffset(offset, cellvefs_dim...)
  dofs_all = IndexCellValueByGlobalAppend(dim_to_nface_eqclass...)
  cell_eqclass = CellVectorByComposition(cellvefs, dofs_all)
  ConformingFESpace{D,Z,T}(
    dim_to_nface_eqclass, cell_eqclass, nfree, nfixed,
    dir_tags, reffe, trian, graph, labels)
end

num_free_dofs(this::ConformingFESpace) = this.num_free_dofs

num_fixed_dofs(this::ConformingFESpace) = this.num_fixed_dofs

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
      if ( _is_fixed(vef_labels[ignf],diri_tags,labels) )
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

@inline function _is_fixed(v,dt,labels)
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
