module ConformingFESpaces

using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery

export ConformingFESpace
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
  _gridgraph::GridGraph
  _labels::FaceLabels
  _basis::CellBasis{Z,T}
end

function ConformingFESpace(
  reffe::LagrangianRefFE{D,T},
  trian::Triangulation{D,Z},
  graph::GridGraph,
  labels::FaceLabels,
  diri_tags::Vector{Int}) where {D,Z,T}
  args = _setup_conforming_fe_fields(reffe,trian,graph,labels,diri_tags,D)
  ConformingFESpace{D,Z,T}(args...)
end

function ConformingFESpace(
  reffe::LagrangianRefFE{D,T},
  trian::Triangulation{D,Z},
  graph::GridGraph,
  labels::FaceLabels) where {D,Z,T}
  return ConformingFESpace(reffe, trian, graph, labels, ())
end

function ConformingFESpace(::Type{T},model::DiscreteModel{D},order,diri_tags) where {D,T}
  grid = Grid(model,D)
  trian = Triangulation(grid)
  graph = GridGraph(model)
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

function interpolate_values(this::ConformingFESpace,f::Function)
  _interpolate_values(this,f)
end

function interpolate_diri_values(this::ConformingFESpace, funs::Vector{<:Function})
  _interpolate_diri_values(this,funs)
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

function _setup_conforming_fe_fields(reffe,trian,graph,labels,diri_tags,D)
  dim_to_nface_eqclass, nfree, ndiri  = _generate_dim_to_nface_to_dofs(
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

function _generate_dim_to_nface_to_dofs(
  reffe::RefFE{D},
  graph::GridGraph,
  labels::FaceLabels,
  diri_tags::Vector{Int}) where D

  i_free_dof = 1
  i_diri_dof = -1

  dim_to_nface_to_dofs = IndexCellVector{Int}[]
  tag_to_labels = labels.tag_to_labels

  for d in 0:D

    nface_to_label = labels_on_dim(labels,d)
    cell_to_nfaces = connections(graph,D,d)
    nface_to_cells = connections(graph,d,D)
    icell = 1
    nface_to_cellowner = get_local_item(nface_to_cells,icell)
    nface_to_lnface = find_local_index(nface_to_cellowner, cell_to_nfaces)

    lnface_to_ldofs = nfacedofs(reffe,d)
    lnface_to_nldofs = length.(lnface_to_ldofs)

    num_nfaces = length(nface_to_cells)

    nface_to_dofs_ptrs = zeros(Int, num_nfaces+1)
    nface_to_dofs_data = Int[]

    i_free_dof, i_diri_dof = _generate_nface_to_dofs!(
      nface_to_dofs_data,
      nface_to_dofs_ptrs,
      nface_to_lnface,
      lnface_to_nldofs,
      nface_to_label,
      diri_tags,
      tag_to_labels,
      i_free_dof,
      i_diri_dof)

    length_to_ptrs!(nface_to_dofs_ptrs)

    nface_to_dofs = CellVectorFromDataAndPtrs(
      nface_to_dofs_data, nface_to_dofs_ptrs)

    push!(dim_to_nface_to_dofs, nface_to_dofs)

  end

  return (dim_to_nface_to_dofs, i_free_dof-1, -i_diri_dof-1)

end

function _generate_nface_to_dofs!(
  nface_to_dofs_data,
  nface_to_dofs_ptrs,
  nface_to_lnface,
  lnface_to_nldofs,
  nface_to_label,
  diri_tags,
  tag_to_labels,
  i_free_dof,
  i_diri_dof)

  for (nface, lnface) in enumerate(nface_to_lnface)

    nldofs = lnface_to_nldofs[lnface]
    label = nface_to_label[nface]
    isdiri = _is_diri(label,diri_tags,tag_to_labels)

    if isdiri
      for i in 1:nldofs
        push!(nface_to_dofs_data,i_diri_dof)
        i_diri_dof += -1
      end
    else
      for i in 1:nldofs
        push!(nface_to_dofs_data,i_free_dof)
        i_free_dof += 1
      end
    end

    nface_to_dofs_ptrs[nface+1] = nldofs
  end

  (i_free_dof, i_diri_dof)
end

#function _generate_dim_to_nface_eqclass(
#  reffe::RefFE{D,T},
#  graph::GridGraph,
#  labels::FaceLabels,
#  diri_tags::Vector{Int}) where {D,T}
#
#  dim_to_nface_eqclass = Int[]
#  c=1
#  c_n = -1
#  for vef_dim in 0:D
#    vefcells= connections(graph,vef_dim,D)
#    cellvefs= connections(graph,D,vef_dim)
#    vef_labels = labels_on_dim(labels,vef_dim)
#    num_vefs = length(vefcells)
#    nfdofs=Array{Array{Int64},1}(undef,num_vefs)
#    nfdofs_l = Int[]
#    nfdofs_g = zeros(Int, num_vefs+1)
#    nfdofs_g[1] = 1
#    for (ignf,nf) in enumerate(vefcells)
#      owner_cell = nf[1]
#      lid_vef_dim = findfirst(i->i==ignf,cellvefs[owner_cell])
#      lid_vef = reffe.polytope.nf_dim[end][vef_dim+1][1]+lid_vef_dim-1
#      # @santiagobadia : Better a method for nfs of a particular type...
#      num_nf_dofs = length(reffe.nfacedofs[lid_vef])
#      if ( _is_diri(vef_labels[ignf],diri_tags,labels.tag_to_labels) )
#        nfdofs_l = Int[nfdofs_l..., c_n:-1:c_n-num_nf_dofs+1... ]
#        c_n -= num_nf_dofs
#      else
#        nfdofs_l = Int[nfdofs_l..., c:c+num_nf_dofs-1... ]
#        c += num_nf_dofs
#      end
#      nfdofs_g[ignf+1] += num_nf_dofs + nfdofs_g[ignf]
#    end
#    dim_to_nface_eqclass = [
#      dim_to_nface_eqclass..., CellVectorFromDataAndPtrs(nfdofs_l, nfdofs_g) ]
#  end
#  return [ dim_to_nface_eqclass , c-1, -c_n-1 ]
#end

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

function _interpolate_values(fesp::ConformingFESpace{D,Z,T},fun::Function) where {D,Z,T}
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

function _interpolate_diri_values(fesp::ConformingFESpace{D,Z,T},funs) where {D,Z,T}
  labels = fesp._labels
  nf_labs_all = [ labels_on_dim(labels,idim) for idim in 0:D]
  nf_dofs_all = fesp.dim_to_nface_eqclass
  dtags = fesp.diri_tags
  @assert length(dtags) == length(funs)
  E = eltype(T)
  diri_dofs_all = zeros(E, num_diri_dofs(fesp))
  for (ifunc,f) in enumerate(funs)
    _ , fh_diri_dofs = interpolate_values(fesp,f)
    # @santiagobadia: For performance issues, implement a new interpolate
    # restricted to cells on the boundary for performance
    for idim in 0:D
      nf_labs = nf_labs_all[idim+1]
      nf_dofs = nf_dofs_all[idim+1]
      # How to extract this part? Do it correctly, with a FEFunction
      for (nf,nflab) in enumerate(nf_labs)
        if (_is_diri(nflab, (dtags[ifunc],), labels.tag_to_labels))
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

@inline function _is_diri(label,diritags,tag_to_labels)
  for tag in diritags
    for dirilabel in tag_to_labels[tag]
      if (label == dirilabel)
        return true
      end
    end
  end
  return false
end

end # module
