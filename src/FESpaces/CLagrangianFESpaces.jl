module CLagrangianFESpaces

# TODO move to TensorValues
using StaticArrays
using TensorValues
mutable(::Type{MultiValue{S,T,N,L}}) where {S,T,N,L} = MArray{S,T,N,L}

using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery
using Gridap.DOFBases: _length
using Gridap.ConformingFESpaces: _CellField

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

export CLagrangianFESpace

struct CLagrangianFESpace{D,Z,T,S,B} <: FESpace{D,Z,T}
  grid::Grid{D,Z}
  node_to_label::Vector{Int}
  tag_to_labels::Vector{Vector{Int}}
  diritags::Vector{Int}
  dirimasks::Vector{B}
  nfreedofs::Int
  ndiridofs::Int
  node_and_comp_to_dof::Vector{S}
  cell_to_dofs::CellVector{Int}
  reffe::LagrangianRefFE{Z,T}
  cellbasis::CellBasis{Z,T}
  dirinodes::Vector{Int}
end

function CLagrangianFESpace(
  ::Type{T},
  grid::Grid{D,Z},
  node_to_label::AbstractVector{<:Integer},
  tag_to_labels::Vector{Vector{Int}},
  diritags::Vector{Int},
  dirimasks::Vector{B}) where {D,Z,T,B}

  node_and_comp_to_dof, nfreedofs, ndiridofs, dirinodes =
    _setup_node_and_comp_to_dof(
    T, node_to_label, tag_to_labels, diritags, dirimasks)

  cell_to_nodes = cells(grid)

  reffe = _setup_reffe(T,celltypes(grid),cellorders(grid))
  lnode_and_comp_to_ldof = reffe.dofbasis.node_and_comp_to_dof

  cell_to_dofs = _setup_cell_to_dofs(
    cell_to_nodes,
    node_and_comp_to_dof,
    lnode_and_comp_to_ldof)

  cellbasis = _setup_cellbasis(reffe,grid)

  S = _S(T)

  CLagrangianFESpace{D,Z,T,S,B}(
    grid,
    node_to_label,
    tag_to_labels,
    diritags,
    dirimasks,
    nfreedofs,
    ndiridofs,
    node_and_comp_to_dof,
    cell_to_dofs,
    reffe,
    cellbasis,
    dirinodes)

end

function CLagrangianFESpace(
  ::Type{T},model::DiscreteModel,order,diritags,dirimasks) where T

  grid = Grid(model)
  _check_order(order,cellorders(grid))
  facelabels = FaceLabels(model)
  node_dim = 0
  node_to_label = labels_on_dim(facelabels,node_dim)
  tag_to_labels = facelabels.tag_to_labels

  CLagrangianFESpace(
    T,grid,node_to_label,tag_to_labels,diritags,dirimasks)

end

num_free_dofs(fesp::CLagrangianFESpace) = fesp.nfreedofs

num_diri_dofs(fesp::CLagrangianFESpace) = fesp.ndiridofs

diri_tags(fesp::CLagrangianFESpace) = fesp.diritags

apply_constraints(
  ::CLagrangianFESpace, cellvec::CellVector, cellids::CellNumber) = cellvec

apply_constraints_rows(
  ::CLagrangianFESpace, cellmat::CellMatrix, cellids::CellNumber) = cellmat

apply_constraints_cols(
  ::CLagrangianFESpace, cellmat::CellMatrix, cellids::CellNumber) = cellmat

celldofids(fesp::CLagrangianFESpace) = fesp.cell_to_dofs

function interpolate_values(fesp::CLagrangianFESpace,fun::Function)
  zh = zero(fesp)
  fdof_to_val = free_dofs(zh)
  ddof_to_val = diri_dofs(zh)
  _fill_interpolated_vals!(fdof_to_val,ddof_to_val,fesp,fun)
  (fdof_to_val, ddof_to_val)
end

function interpolate_diri_values(fesp::CLagrangianFESpace, funs::Vector{<:Function})
  
  T = value_type(fesp)
  E = eltype(T)
  ndiri = num_diri_dofs(fesp)
  ddof_to_val = zeros(E,ndiri)
  node_to_coords = points(fesp.grid)
  node_and_comp_to_dof = fesp.node_and_comp_to_dof
  node_to_label = fesp.node_to_label
  tag_to_labels = fesp.tag_to_labels
  dirinodes = fesp.dirinodes

  for (fun,diritag,dirimask) in zip(funs,fesp.diritags,fesp.dirimasks)
    _fill_diri_values!(
      ddof_to_val,
      fun,
      diritag,
      dirimask,
      dirinodes,
      node_to_label,
      tag_to_labels,
      node_to_coords,
      node_and_comp_to_dof)
  end

  ddof_to_val
end

function CellField(
  fesp::CLagrangianFESpace{D,Z,T},
  free_dofs::AbstractVector{E},
  diri_dofs::AbstractVector{E}) where {D,Z,T,E}

  _CellField(fesp, free_dofs, diri_dofs, T, E)
end

CellBasis(fesp::CLagrangianFESpace{D,Z,T}) where {D,Z,T} = fesp.cellbasis

# Helpers


function _setup_node_and_comp_to_dof(
  ::Type{T}, node_to_label, tag_to_labels, diritags, dirimasks) where T

  S = _S(T)
  nnodes = length(node_to_label)
  node_and_comp_to_dof = zeros(S,nnodes)
  dirinodes = Int[]

  nf, nd = _fill_node_and_comp_to_dof!(
    node_and_comp_to_dof, dirinodes, node_to_label, tag_to_labels, diritags, dirimasks)

  (node_and_comp_to_dof, nf, nd, dirinodes)

end

function _fill_node_and_comp_to_dof!(
  node_and_comp_to_dof::Vector{S},
  dirinodes,
  node_to_label,
  tag_to_labels,
  diritags,
  dirimasks) where S

  nnodes = length(node_to_label)
  fdof = 1
  ddof = -1

  for node in 1:nnodes
    label = node_to_label[node]

    comp_to_dof, fdof, ddof, isdirinode = 
      _compute_comp_to_dof(S,label,tag_to_labels,diritags,dirimasks,fdof,ddof)

    if isdirinode
      push!(dirinodes,node)
    end

    node_and_comp_to_dof[node] = comp_to_dof

  end

  (fdof-1,-(ddof+1))

end

function _compute_comp_to_dof(
  ::Type{S}, label, tag_to_labels, diritags, dirimasks, fdof, ddof) where S <: MultiValue

  comp_to_dof = zero(mutable(S))
  comps = eachindex(comp_to_dof)

  isdirinode = false
  for comp in comps
    dof, fdof, ddof = _compute_dof(label,tag_to_labels,diritags,dirimasks,comp,fdof,ddof)
    comp_to_dof[comp] = dof
    if dof < 0
      isdirinode = true
    end
  end

  (comp_to_dof, fdof, ddof, isdirinode)

end

function _compute_comp_to_dof(
  ::Type{S}, label, tag_to_labels, diritags, dirimasks, fdof, ddof) where S <: Real

  comp = 1
  dof, fdof, ddof =  _compute_dof(
    label,tag_to_labels,diritags,dirimasks,comp,fdof,ddof)

  isdirinode = dof < 0

  (dof, fdof, ddof, isdirinode)

end

function _compute_dof(label,tag_to_labels,diritags,dirimasks,comp,fdof,ddof)

  isfree = true
  for (i,diritag) in enumerate(diritags)

    if label in tag_to_labels[diritag]
      comp_to_mask = dirimasks[i]
      mask = comp_to_mask[comp]
      if mask
        dof = ddof
        ddof -= 1
        isfree = false
        break
      end
    end
  end

  if isfree
    dof = fdof
    fdof += 1
  end

  (dof,fdof,ddof)

end


_S(::Type{<:Real}) = Int

_S(::Type{MultiValue{S,T,N,L}}) where {S,T,N,L} = MultiValue{S,Int,N,L}

function _setup_reffe(ct,co)
  @notimplemented
end

function _setup_reffe(
  ::Type{T},
  ct::ConstantCellValue{<:NTuple{Z}},
  co::ConstantCellValue) where {T,Z}
  t = ct.value
  o = co.value
  p = Polytope(t)
  LagrangianRefFE{Z,T}(p,fill(o,Z))
end

function _setup_cellbasis(reffe,grid)
  lshapes = shfbasis(reffe)
  nc = ncells(grid)
  shb = ConstantCellValue(lshapes, nc)
  trian = Triangulation(grid)
  phi = CellGeomap(trian)
  basis = attachgeomap(shb,phi)
  basis
end

function _setup_cell_to_dofs(
  cell_to_nodes,
  node_and_comp_to_dof::Vector{S},
  lnode_and_comp_to_ldof) where S

  ncomp = _length(S)
  ncells = length(cell_to_nodes)
  nlnodes, ncomp = size(lnode_and_comp_to_ldof)
  nldofs = nlnodes * ncomp

  ndata = ncells*nldofs
  I = eltype(S)

  cell_to_dofs_data = zeros(I,ndata)

  _fill_cell_to_dofs!(
    cell_to_dofs_data,
    cell_to_nodes,
    node_and_comp_to_dof,
    lnode_and_comp_to_ldof)

  CellVectorFromDataAndStride(cell_to_dofs_data,nldofs)

end

function  _fill_cell_to_dofs!(
  cell_to_dofs_data,
  cell_to_nodes,
  node_and_comp_to_dof,
  lnode_and_comp_to_ldof)

  k = 0
  
  nlnodes, ncomps = size(lnode_and_comp_to_ldof)
  nldofs = nlnodes * ncomps

  for nodes in cell_to_nodes

    for (lnode,node) in enumerate(nodes)
      for comp in 1:ncomps
        ldof = lnode_and_comp_to_ldof[lnode,comp]
        dof = node_and_comp_to_dof[node][comp]
        cell_to_dofs_data[k+ldof] = dof
      end
    end

    k += nldofs

  end

end

function _fill_interpolated_vals!(fdof_to_val,ddof_to_val,fesp,fun)

  x = points(fesp.grid)
  fx = fun.(x)
  T = value_type(fesp)
  @assert eltype(fx) == T
  z = zero(T)
  comps = eachindex(z)

  _fill_interpolated_vals_kernel!(
    fdof_to_val,
    ddof_to_val,
    fesp.node_and_comp_to_dof,
    fx,
    comps)

end

function  _fill_interpolated_vals_kernel!(
  fdof_to_val,
  ddof_to_val,
  node_and_comp_to_dof,
  node_and_comp_to_val,
  comps)

  nnodes = length(node_and_comp_to_dof)
  @assert nnodes == length(node_and_comp_to_val)

  for node in 1:nnodes
    comp_to_dof = node_and_comp_to_dof[node]
    comp_to_val = node_and_comp_to_val[node]
    for comp in comps
      dof = comp_to_dof[comp]
      val = comp_to_val[comp]
      if dof > 0
        fdof_to_val[dof] = val
      elseif dof < 0
        ddof_to_val[-dof] = val
      else
        @unreachable
      end
    end
  end

end

function    _fill_diri_values!(
  ddof_to_val,
  fun,
  diritag,
  dirimask,
  dirinodes,
  node_to_label,
  tag_to_labels,
  node_to_coords,
  node_and_comp_to_dof)

  for node in dirinodes

    label = node_to_label[node]
    if label in tag_to_labels[diritag]
      x = node_to_coords[node]
      fx = fun(x)
      comp_to_dof = node_and_comp_to_dof[node]
      for comp in eachindex(fx)
        mask = dirimask[comp]
        if mask
          val = fx[comp]
          dof = comp_to_dof[comp]
          ddof = -dof
          ddof_to_val[ddof] = val
        end
      end
    end

  end

end

function _check_order(order,co)
  @notimplemented
end

function _check_order(order,co::ConstantCellValue)
  @notimplementedif order != co.value
end

end # module
