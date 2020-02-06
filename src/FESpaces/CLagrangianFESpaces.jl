
struct CLagrangianFESpace{S} <: SingleFieldFESpace
  space::UnsconstrainedFESpace
  grid::Grid
  dof_to_node::Vector{Int}
  dof_to_comp::Vector{Int8}
  node_and_comp_to_dof::Vector{S}

  function CLagrangianFESpace(::Type{T},grid::Grid) where T
    space, dof_to_node, dof_to_comp, node_and_comp_to_dof = _generate_clargangian_fespace(T,grid)
    S = eltype(node_and_comp_to_dof)
    new{S}(space,grid,dof_to_node,dof_to_comp,node_and_comp_to_dof)
  end
end

# FESpace interface

function num_free_dofs(f::CLagrangianFESpace)
  num_free_dofs(f.space)
end

function get_cell_basis(f::CLagrangianFESpace)
  get_cell_basis(f.space)
end

function zero_free_values(::Type{T},f::CLagrangianFESpace) where T
  zero_free_values(T,f.space)
end

function apply_constraints_matrix_cols(f::CLagrangianFESpace,cellmat,cellids)
  apply_constraints_matrix_cols(f.space,cellmat,cellids)
end

function apply_constraints_matrix_rows(f::CLagrangianFESpace,cellmat,cellids)
  apply_constraints_matrix_rows(f.space,cellmat,cellids)
end

function apply_constraints_vector(f::CLagrangianFESpace,cellvec,cellids)
  apply_constraints_vector(f.space,cellvec,cellids)
end

# SingleFieldFESpace interface

function get_cell_dofs(f::CLagrangianFESpace)
  get_cell_dofs(f.space)
end

function get_cell_dof_basis(f::CLagrangianFESpace)
  get_cell_dof_basis(f.space)
end

function num_dirichlet_dofs(f::CLagrangianFESpace)
  num_dirichlet_dofs(f.space)
end

function num_dirichlet_tags(f::CLagrangianFESpace)
  num_dirichlet_tags(f.space)
end

function zero_dirichlet_values(f::CLagrangianFESpace)
  zero_dirichlet_values(f.space)
end

function get_dirichlet_dof_tag(f::CLagrangianFESpace)
  get_dirichlet_dof_tag(f.space)
end

function scatter_free_and_dirichlet_values(f::CLagrangianFESpace,free_values,dirichlet_values)
  scatter_free_and_dirichlet_values(f.space,free_values,dirichlet_values)
end

function gather_free_and_dirichlet_values(f::CLagrangianFESpace,cell_vals)
  gather_free_and_dirichlet_values(f.space,cell_vals)
end

function gather_dirichlet_values(f::CLagrangianFESpace,cell_vals)
  gather_dirichlet_values(f.space,cell_vals)
end

# Helpers

function _generate_clargangian_fespace(T,grid)

    grid_reffes = get_reffes(grid)
    function fun(reffe)
      p = get_polytope(reffe)
      orders = get_orders(reffe)
      LagrangianRefFE(T,p,orders)
    end
    reffes = map(fun,grid_reffes)

    nnodes = num_nodes(grid)
    dof_to_node, dof_to_comp, node_and_comp_to_dof = _generate_dof_layout_component_major(T,nnodes)

    cell_nodes = get_cell_nodes(grid)
    cell_to_ctype = get_cell_type(grid)

    cell_dofs = _generate_cell_dofs_clagrangian_fespace(
      T,
      cell_nodes,
      reffes,
      cell_to_ctype,
      node_and_comp_to_dof)

    cell_map = get_cell_map(grid)
    cell_shapefuns, cell_dof_basis = compute_cell_space(reffes, cell_to_ctype, cell_map)

    nfree = length(dof_to_node)
    ndirichlet = 0

    dirichlet_dof_tag = Int8[]
    dirichlet_cells = Int[]
    ntags = 0

    space = UnsconstrainedFESpace(
      nfree,
      ndirichlet,
      cell_dofs,
      cell_shapefuns,
      cell_dof_basis,
      cell_map,
      dirichlet_dof_tag,
      dirichlet_cells,
      ntags)

    space, dof_to_node, dof_to_comp, node_and_comp_to_dof

end

function _generate_dof_layout_component_major(::Type{<:Real},nnodes::Integer)
  ndofs = nnodes
  dof_to_comp = ones(Int8,ndofs)
  dof_to_node = collect(1:nnodes)
  node_and_comp_to_dof = dof_to_node
  (dof_to_node, dof_to_comp, node_and_comp_to_dof)
end

function _generate_dof_layout_component_major(::Type{T},nnodes::Integer) where T
  ncomps = n_components(T)
  V = change_eltype(T,Int)
  ndofs = ncomps*nnodes
  dof_to_comp = zeros(Int8,ndofs)
  dof_to_node = zeros(Int,ndofs)
  node_and_comp_to_dof = zeros(V,nnodes)
  m = zero(mutable(V))
  for node in 1:nnodes
    o = ncomps*(node-1)
    for comp in 1:ncomps
      dof = comp+o
      dof_to_comp[dof] = comp
      dof_to_node[dof] = node
      m[comp] = dof
    end
    node_and_comp_to_dof[node] = m
  end
  (dof_to_node, dof_to_comp, node_and_comp_to_dof)
end

function _generate_cell_dofs_clagrangian_fespace(
  ::Type{<:Real},
  cell_nodes,
  reffes,
  cell_to_ctype,
  node_and_comp_to_dof)

  cell_nodes
end

function _generate_cell_dofs_clagrangian_fespace(
  ::Type{T},
  cell_nodes,
  reffes,
  cell_to_ctype,
  node_and_comp_to_dof) where T

  ncomps = n_components(T)

  ctype_to_lnode_to_comp_to_ldof = map(get_node_and_comp_to_dof,reffes)
  ctype_to_num_ldofs = map(num_dofs,reffes)

  cell_dofs = _allocate_cell_dofs_clagrangian_fespace(ctype_to_num_ldofs,cell_to_ctype)

  cache = array_cache(cell_nodes)
  _fill_cell_dofs_clagrangian_fespace!(
    cell_dofs,
    cache,
    cell_nodes,
    ctype_to_lnode_to_comp_to_ldof,
    cell_to_ctype,
    node_and_comp_to_dof,
    ncomps)

  cell_dofs
end

function _allocate_cell_dofs_clagrangian_fespace(ctype_to_num_ldofs,cell_to_ctype)
  ncells = length(cell_to_ctype)
  ptrs = zeros(Int32,ncells+1)
  for cell in 1:ncells
    ctype = cell_to_ctype[cell]
    nldofs = ctype_to_num_ldofs[ctype]
    ptrs[cell+1] = nldofs
  end
  length_to_ptrs!(ptrs)
  ndata = ptrs[end]-1
  data = zeros(Int,ndata)
  Table(data,ptrs)
end

function _fill_cell_dofs_clagrangian_fespace!(
  cell_dofs,
  cache,
  cell_nodes,
  ctype_to_lnode_and_comp_to_ldof,
  cell_to_ctype,
  node_and_comp_to_dof,
  ncomps)

  ncells = length(cell_to_ctype)

  for cell in 1:ncells
    nodes = getindex!(cache,cell_nodes,cell)
    ctype = cell_to_ctype[cell]
    lnode_and_comp_to_ldof = ctype_to_lnode_and_comp_to_ldof[ctype]
    p = cell_dofs.ptrs[cell]-1
    for (lnode, node) in enumerate(nodes)
      for comp in 1:ncomps
        ldof = lnode_and_comp_to_ldof[lnode][comp]
        dof = node_and_comp_to_dof[node][comp]
        cell_dofs.data[p+ldof] = dof
      end
    end
  end

end

