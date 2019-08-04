module CLagrangianFESpaces

using Gridap
using Gridap.DOFBases: _length
using Gridap.ConformingFESpaces: _CellField

struct CLagrangianFESpace{D,Z,T,S,B} <: FESpace{D,Z,T}
  grid::Grid{D,Z}
  node_to_tag::Vector{Int}
  diritags::Vector{Int}
  dirimasks::Vector{B}
  nfreedofs::Int
  ndiridofs::Int
  node_and_comp_to_dof::Vector{S}
  cell_to_dofs::CellVector{Int}
  reffe::LagrangianRefFE{Z,T}
  cellbasis::CellBasis{Z,T}
end

function CLagrangianFESpace(
  ::Type{T},
  grid::Grid{D,Z},
  node_to_tag::Vector{Int},
  diritags::Vector{Int},
  dirimasks::Vector{B}) where {D,Z,T,B}

  node_and_comp_to_dof = _setup_node_and_comp_to_dof(
    node_to_tag, diritags, dirimasks)

  cell_to_nodes = cells(grid)

  reffe = _setup_reffe(celltypes(grid),cellorders(grid))
  lnode_and_comp_to_ldof = reffe.dofbasis.node_and_comp_to_dof

  cell_to_dofs = _setup_cell_to_dofs(
    cell_to_nodes,
    node_and_comp_to_dof,
    lnode_and_comp_to_ldof)

  cellbasis = _setup_cellbasis(reffe,grid)






  p = polytope(grid_reffe)

  reffe = LagrangianRefFE{Z,T}(p,fill(order,Z))

  lshapes = shfbasis(reffe)

  nc = length(cell_to_nodes)
  shb = ConstantCellValue(lshapes, nc)
  trian = Triangulation(grid)
  phi = CellGeomap(trian)
  basis = attachgeomap(shb,phi)

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
  zh = zero(fesp)
  fdof_to_val = free_dofs(zh)
  ddof_to_val = diri_dofs(zh)
  for fun in funs
    _fill_interpolated_vals!(fdof_to_val,ddof_to_val,fesp,fun)
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
  ::Type{T}, node_to_tag, diritags, dirimasks) where T
  S = _S(T)
  nnodes = length(node_to_tag)
  node_and_comp_to_dof = zeros(S,nnodes)
  comp_to_dof = zero(mutable(S))
  z = zero(eltype(S))
  _fill_node_and_comp_to_dof!(
    node_and_comp_to_dof, node_to_tag, diritags, dirimasks, comp_to_dof, z)
  node_and_comp_to_dof
end

function _fill_node_and_comp_to_dof!(
  node_and_comp_to_dof, node_to_tag, diritags, dirimasks, comp_to_dof, z)

  nnodes = length(node_to_tag)
  comps = eachindex(comp_to_dof)
  fdof = 1
  ddof = -1

  for node in 1:nnodes
    tag = node_to_tag[node]

    for comp in comps
      comp_to_dof[comp] = z


      isfree = true
      for (i,diritag) in enumerate(diritags)
        if tag == diritag
          comp_to_mask = dirimasks[i]
          mask = comp_to_mask[comp]
          if mask
            comp_to_dof[comp] = ddof
            ddof -= 1
            isfree = false
            break
          end
        end
      end

      if isfree
        comp_to_dof[comp] = fdof
        fdof += 1
      end

      node_and_comp_to_dof[node] = comp_to_dof

    end

  end

end

_S(::Type{<:Real}) = Int

_S(::Type{MultiValue{S,T,N,L}}) where {S,T,N,L} = MultiValue{S,Int,N,L}

# TODO move to TensorValues
mutable(::Type{MultiValue{S,T,N,L}}) where {S,T,N,L} = MArray{S,T,N,L}

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
        fdof_to_val[-dof] = val
      else
        @unreachable
      end
    end
  end

end

end # module
