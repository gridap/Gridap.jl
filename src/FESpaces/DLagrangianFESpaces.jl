module DLagrangianFESpaces

using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery
using Gridap.CLagrangianFESpaces: _setup_reffe
using Gridap.CLagrangianFESpaces: _setup_cellbasis
using Gridap.CLagrangianFESpaces: _S
using Gridap.CLagrangianFESpaces: _compute_comp_to_dof
using Gridap.CLagrangianFESpaces: _setup_grid
using Gridap.ConformingFESpaces: _CellField

export DLagrangianFESpace

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

struct DLagrangianFESpace{D,Z,T,B} <: FESpace{D,Z,T}
  grid::Grid{D,Z}
  node_to_label::Vector{Int}
  tag_to_labels::Vector{Vector{Int}}
  diritags::Vector{Int}
  dirimasks::Vector{B}
  nfreedofs::Int
  ndiridofs::Int
  cell_to_dofs::CellVector{Int}
  reffe::LagrangianRefFE{Z,T}
  cellbasis::CellBasis{Z,T}
  diricells::Vector{Int}
end

function DLagrangianFESpace(
  ::Type{T},
  grid::Grid{D,Z},
  node_to_label::AbstractVector{<:Integer},
  tag_to_labels::Vector{Vector{Int}},
  diritags::Vector{Int},
  dirimasks::Vector{B}) where {D,Z,T,B}

  cell_to_nodes = cells(grid)
  reffe = _setup_reffe(T,celltypes(grid),cellorders(grid))
  cellbasis = _setup_cellbasis(reffe,grid)

  lnode_and_comp_to_ldof = reffe.dofbasis.node_and_comp_to_dof

  cell_to_dofs, nfreedofs, ndiridofs, diricells = _setup_cell_to_dofs(
    cell_to_nodes,
    node_to_label,
    tag_to_labels,
    diritags,
    dirimasks,
    lnode_and_comp_to_ldof,
    T)

  DLagrangianFESpace{D,Z,T,B}(
    grid,
    node_to_label,
    tag_to_labels,
    diritags,
    dirimasks,
    nfreedofs,
    ndiridofs,
    cell_to_dofs,
    reffe,
    cellbasis,
    diricells)

end

function DLagrangianFESpace(
  ::Type{T},model::DiscreteModel,order,diritags,dirimasks) where T

  grid, node_to_label, tag_to_labels = _setup_grid(model,order)

  DLagrangianFESpace(
    T,grid,node_to_label,tag_to_labels,diritags,dirimasks)

end

num_free_dofs(fesp::DLagrangianFESpace) = fesp.nfreedofs

num_diri_dofs(fesp::DLagrangianFESpace) = fesp.ndiridofs

diri_tags(fesp::DLagrangianFESpace) = fesp.diritags

apply_constraints(
  ::DLagrangianFESpace, cellvec::CellVector, cellids::CellNumber) = cellvec

apply_constraints_rows(
  ::DLagrangianFESpace, cellmat::CellMatrix, cellids::CellNumber) = cellmat

apply_constraints_cols(
  ::DLagrangianFESpace, cellmat::CellMatrix, cellids::CellNumber) = cellmat

celldofids(fesp::DLagrangianFESpace) = fesp.cell_to_dofs

function CellField(
  fesp::DLagrangianFESpace{D,Z,T},
  free_dofs::AbstractVector{E},
  diri_dofs::AbstractVector{E}) where {D,Z,T,E}

  _CellField(fesp, free_dofs, diri_dofs, T, E)
end

function interpolate_values(fesp::DLagrangianFESpace,fun::Function)

  zh = zero(fesp)
  fdof_to_val = free_dofs(zh)
  ddof_to_val = diri_dofs(zh)
  node_to_coords = points(fesp.grid)
  lnode_and_comp_to_ldof = fesp.reffe.dofbasis.node_and_comp_to_dof
  cell_to_dofs = fesp.cell_to_dofs
  cell_to_nodes = cells(fesp.grid)
  node_comp_to_val = fun.(node_to_coords)

  _fill_interpolated_vals!(
    fdof_to_val,
    ddof_to_val,
    node_comp_to_val,
    cell_to_nodes,
    cell_to_dofs,
    lnode_and_comp_to_ldof)

  (fdof_to_val, ddof_to_val)
end

function interpolate_diri_values(fesp::DLagrangianFESpace, funs::Vector{<:Function})

  T = value_type(fesp)
  E = eltype(T)
  ndiri = num_diri_dofs(fesp)
  ddof_to_val = zeros(E,ndiri)
  node_to_label = fesp.node_to_label
  tag_to_labels = fesp.tag_to_labels
  diricells = fesp.diricells
  lnode_and_comp_to_ldof = fesp.reffe.dofbasis.node_and_comp_to_dof
  cell_to_dofs = fesp.cell_to_dofs
  cell_to_nodes = cells(fesp.grid)
  node_to_coords = points(fesp.grid)

  for (fun,diritag,dirimask) in zip(funs,fesp.diritags,fesp.dirimasks)
    _fill_diri_values!(
      ddof_to_val,
      fun,
      diritag,
      dirimask,
      diricells,
      node_to_label,
      tag_to_labels,
      node_to_coords,
      cell_to_nodes,
      cell_to_dofs,
      lnode_and_comp_to_ldof)
  end

  ddof_to_val

end

CellBasis(fesp::DLagrangianFESpace{D,Z,T}) where {D,Z,T} = fesp.cellbasis

Triangulation(fesp::DLagrangianFESpace) = Triangulation(fesp.grid)

# Helpers

function _fill_diri_values!(
  ddof_to_val,
  fun,
  diritag,
  dirimask,
  diricells,
  node_to_label,
  tag_to_labels,
  node_to_coords,
  cell_to_nodes,
  cell_to_dofs,
  lnode_and_comp_to_ldof)

  for cell in diricells

    lnode_to_node = cell_to_nodes[cell]
    ldof_to_dof = cell_to_dofs[cell]


    for (lnode,node) in enumerate(lnode_to_node)
      label = node_to_label[node]
      if label in tag_to_labels[diritag]

        coords = node_to_coords[node]
        comp_to_val = fun(coords)

        for comp in eachindex(dirimask)
          mask = dirimask[comp]
          if mask
            val = comp_to_val[comp]
            ldof = lnode_and_comp_to_ldof[lnode,comp]
            dof = ldof_to_dof[ldof]
            ddof = -dof
            ddof_to_val[ddof] = val
          end
        end

      end
    end

  end

end

function  _fill_interpolated_vals!(
  fdof_to_val,
  ddof_to_val,
  node_comp_to_val,
  cell_to_nodes,
  cell_to_dofs,
  lnode_and_comp_to_ldof)

  ncells = length(cell_to_nodes)

  nlnodes, ncomps = size(lnode_and_comp_to_ldof)

  for cell in 1:ncells

    lnode_to_node = cell_to_nodes[cell]
    ldof_to_dof = cell_to_dofs[cell]

    for (lnode,node) in enumerate(lnode_to_node)
      comp_to_val = node_comp_to_val[node]
      for comp in 1:ncomps
        ldof = lnode_and_comp_to_ldof[lnode,comp]
        dof = ldof_to_dof[ldof]
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

end

function _setup_cell_to_dofs(
  cell_to_nodes,
  node_to_label,
  tag_to_labels,
  diritags,
  dirimasks,
  lnode_and_comp_to_ldof,
  ::Type{T}) where T

  ncells = length(cell_to_nodes)
  nlnodes, ncomp = size(lnode_and_comp_to_ldof)
  nldofs = nlnodes * ncomp

  ndata = ncells*nldofs

  cell_to_dofs_data = zeros(Int,ndata)
  diricells = Int[]

  nfreedofs, ndiridofs = _fill_cell_to_dofs!(
    cell_to_dofs_data,
    diricells,
    cell_to_nodes,
    node_to_label,
    tag_to_labels,
    diritags,
    dirimasks,
    lnode_and_comp_to_ldof,
    T)

  cell_to_dofs = CellVectorFromDataAndStride(cell_to_dofs_data,nldofs)

  (cell_to_dofs, nfreedofs, ndiridofs, diricells)
end

function _fill_cell_to_dofs!(
  cell_to_dofs_data,
  diricells,
  cell_to_nodes,
  node_to_label,
  tag_to_labels,
  diritags,
  dirimasks,
  lnode_and_comp_to_ldof,
  ::Type{T}) where T

  ncells = length(cell_to_nodes)
  nlnodes, ncomp = size(lnode_and_comp_to_ldof)
  nldofs = nlnodes * ncomp
  S = _S(T)

  fdof = 1
  ddof = -1
  k = 0

  for cell in 1:ncells

    nodes = cell_to_nodes[cell]
    isdiricell = false

    for (lnode,node) in enumerate(nodes)

      label = node_to_label[node]
      comp_to_dof, fdof, ddof, isdirinode = 
        _compute_comp_to_dof(
        S,label,tag_to_labels,diritags,dirimasks,fdof,ddof)

      for comp in eachindex(comp_to_dof)
        ldof = lnode_and_comp_to_ldof[lnode,comp]
        dof = comp_to_dof[comp]
        cell_to_dofs_data[k+ldof] = dof
      end

      if isdirinode
        isdiricell = true
      end

    end

    if isdiricell
      push!(diricells,cell)
    end

    k += nldofs

  end

  (fdof-1, -(ddof+1))

end




end # module
