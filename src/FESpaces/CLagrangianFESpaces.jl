
"""
    struct NodeToDofGlue{T}
      free_dof_to_node::Vector{Int32}
      free_dof_to_comp::Vector{Int16}
      dirichlet_dof_to_node::Vector{Int32}
      dirichlet_dof_to_comp::Vector{Int16}
      node_and_comp_to_dof::Vector{T}
    end
"""
struct NodeToDofGlue{T}
  free_dof_to_node::Vector{Int32}
  free_dof_to_comp::Vector{Int16}
  dirichlet_dof_to_node::Vector{Int32}
  dirichlet_dof_to_comp::Vector{Int16}
  node_and_comp_to_dof::Vector{T}
end

"""
    struct CLagrangianFESpace{T} <: SingleFieldFESpace
      grid::Grid
      glue::NodeToDofGlue{T}
      # + private fields
    end
"""
struct CLagrangianFESpace{S} <: SingleFieldFESpace
  space::UnconstrainedFESpace
  grid::Grid
  glue::NodeToDofGlue{S}
end

"""
    CLagrangianFESpace(::Type{T},grid::Grid) where T
"""
function CLagrangianFESpace(::Type{T},grid::Grid) where T
  vector_type = Vector{_dof_type(T)}
  node_to_tag = fill(Int8(UNSET),num_nodes(grid))
  tag_to_mask = fill(_default_mask(T),0)
  CLagrangianFESpace(T,grid,vector_type,node_to_tag,tag_to_mask)
end

"""
    CLagrangianFESpace(
    ::Type{T},
    grid::Grid,
    vector_type::Type,
    node_to_tag::AbstractVector,
    tag_to_mask::AbstractVector) where T
"""
function CLagrangianFESpace(
  ::Type{T},
  grid::Grid,
  vector_type::Type,
  node_to_tag::AbstractVector{<:Integer},
  tag_to_masks::AbstractVector) where T

  z = zero(T)
  glue, dirichlet_dof_tag = _generate_node_to_dof_glue_component_major(
    z,node_to_tag,tag_to_masks)
  cell_dofs_ids, cell_reffe = _generate_cell_data_clagrangian(z,grid,glue)
  cell_shapefuns, cell_dof_basis = compute_cell_space(cell_reffe, grid)
  cell_is_dirichlet = collect(lazy_map(dofs->any(dof->dof<0,dofs),cell_dofs_ids))

  nfree = length(glue.free_dof_to_node)
  ndirichlet = length(glue.dirichlet_dof_to_node)
  ntags = length(tag_to_masks)
  dirichlet_cells = collect(Int32,findall(cell_is_dirichlet))

  space = UnconstrainedFESpace(
    vector_type,
    nfree,
    ndirichlet,
    cell_dofs_ids,
    cell_shapefuns,
    cell_dof_basis,
    cell_is_dirichlet,
    dirichlet_dof_tag,
    dirichlet_cells,
    ntags)

  CLagrangianFESpace(space,grid,glue)
end

# FESpace interface

ConstraintStyle(::Type{<:CLagrangianFESpace}) = UnConstrained()

function get_free_dof_ids(f::CLagrangianFESpace)
  get_free_dof_ids(f.space)
end

function get_cell_shapefuns(f::CLagrangianFESpace)
  get_cell_shapefuns(f.space)
end

function get_cell_shapefuns_trial(f::CLagrangianFESpace)
  get_cell_shapefuns_trial(f.space)
end

get_vector_type(f::CLagrangianFESpace) = get_vector_type(f.space)

get_dof_value_type(f::CLagrangianFESpace) = get_dof_value_type(f.space)

get_triangulation(f::CLagrangianFESpace) = get_triangulation(f.space)

# SingleFieldFESpace interface

function get_cell_dof_ids(f::CLagrangianFESpace)
  get_cell_dof_ids(f.space)
end

function get_cell_dof_basis(f::CLagrangianFESpace)
  get_cell_dof_basis(f.space)
end

function get_dirichlet_dof_ids(f::CLagrangianFESpace)
  get_dirichlet_dof_ids(f.space)
end

function num_dirichlet_tags(f::CLagrangianFESpace)
  num_dirichlet_tags(f.space)
end

function get_dirichlet_dof_tag(f::CLagrangianFESpace)
  get_dirichlet_dof_tag(f.space)
end

function scatter_free_and_dirichlet_values(f::CLagrangianFESpace,free_values,dirichlet_values)
  scatter_free_and_dirichlet_values(f.space,free_values,dirichlet_values)
end

function gather_free_and_dirichlet_values!(free_values,dirichlet_values,f::CLagrangianFESpace,cell_vals)
  gather_free_and_dirichlet_values!(free_values,dirichlet_values,f.space,cell_vals)
end

function gather_dirichlet_values!(dirichlet_vals,f::CLagrangianFESpace,cell_vals)
    gather_dirichlet_values!(dirichlet_vals,f.space,cell_vals)
end

# Helpers

_default_mask(::Type) = true
_default_mask(::Type{T}) where T <: MultiValue = ntuple(i->true,Val{length(T)}())

_dof_type(::Type{T}) where T = T
_dof_type(::Type{T}) where T<:MultiValue = eltype(T)

function _generate_node_to_dof_glue_component_major(z,node_to_tag,tag_to_mask)
  @notimplementedif length(z) != 1
  nfree_dofs = 0
  ndiri_dofs = 0
  for tag in node_to_tag
    if tag == UNSET
      nfree_dofs += 1
    else
      mask = tag_to_mask[tag]
      if mask
        ndiri_dofs += 1
      else
        nfree_dofs += 1
      end
    end
  end
  free_dof_to_node = Vector{Int32}(undef,nfree_dofs)
  free_dof_to_comp = ones(Int16,nfree_dofs)
  diri_dof_to_node = Vector{Int32}(undef,ndiri_dofs)
  diri_dof_to_comp = ones(Int16,ndiri_dofs)
  diri_dof_to_tag = ones(Int8,ndiri_dofs)
  nnodes = length(node_to_tag)
  node_and_comp_to_dof = Vector{Int32}(undef,nnodes)
  nfree_dofs = 0
  ndiri_dofs = 0
  for (node,tag) in enumerate(node_to_tag)
    if tag == UNSET
      nfree_dofs += 1
      free_dof_to_node[nfree_dofs] = node
      node_and_comp_to_dof[node] = nfree_dofs
    else
      mask = tag_to_mask[tag]
      if mask
        ndiri_dofs += 1
        diri_dof_to_node[ndiri_dofs] = node
        diri_dof_to_tag[ndiri_dofs] = tag
        node_and_comp_to_dof[node] = -ndiri_dofs
      else
        nfree_dofs += 1
        free_dof_to_node[nfree_dofs] = node
        node_and_comp_to_dof[node] = nfree_dofs
      end
    end
  end
  glue = NodeToDofGlue(
   free_dof_to_node,
   free_dof_to_comp,
   diri_dof_to_node,
   diri_dof_to_comp,
   node_and_comp_to_dof)

  glue, diri_dof_to_tag
end

function _generate_node_to_dof_glue_component_major(
  z::MultiValue,node_to_tag,tag_to_masks)
  nfree_dofs = 0
  ndiri_dofs = 0
  ncomps = length(z)
  @assert length(testitem(tag_to_masks)) == ncomps
  for (node,tag) in enumerate(node_to_tag)
    if tag == UNSET
      nfree_dofs += ncomps
    else
      for mask in tag_to_masks[tag]
        if mask
          ndiri_dofs += 1
        else
          nfree_dofs += 1
        end
      end
    end
  end
  free_dof_to_node = Vector{Int32}(undef,nfree_dofs)
  free_dof_to_comp = ones(Int16,nfree_dofs)
  diri_dof_to_node = Vector{Int32}(undef,ndiri_dofs)
  diri_dof_to_comp = ones(Int16,ndiri_dofs)
  diri_dof_to_tag = ones(Int8,ndiri_dofs)
  T = change_eltype(z,Int32)
  nnodes = length(node_to_tag)
  node_and_comp_to_dof = zeros(T,nnodes)
  nfree_dofs = 0
  ndiri_dofs = 0
  m = zero(Mutable(T))
  for (node,tag) in enumerate(node_to_tag)
    if tag == UNSET
      for comp in 1:ncomps
        nfree_dofs += 1
        free_dof_to_node[nfree_dofs] = node
        free_dof_to_comp[nfree_dofs] = comp
        m[comp] = nfree_dofs
      end
    else
      masks = tag_to_masks[tag]
      for comp in 1:ncomps
        mask = masks[comp]
        if mask
          ndiri_dofs += 1
          diri_dof_to_node[ndiri_dofs] = node
          diri_dof_to_comp[ndiri_dofs] = comp
          diri_dof_to_tag[ndiri_dofs] = tag
          m[comp] = -ndiri_dofs
        else
          nfree_dofs += 1
          free_dof_to_node[nfree_dofs] = node
          free_dof_to_comp[nfree_dofs] = comp
          m[comp] = nfree_dofs
        end
      end
    end
    node_and_comp_to_dof[node] = m
  end
  glue = NodeToDofGlue(
   free_dof_to_node,
   free_dof_to_comp,
   diri_dof_to_node,
   diri_dof_to_comp,
   node_and_comp_to_dof)

  glue, diri_dof_to_tag
end

function _generate_cell_data_clagrangian(z,grid,glue)
  ctype_to_reffe = get_reffes(grid)
  cell_to_ctype = get_cell_type(grid)
  function fun(reffe)
    p = get_polytope(reffe)
    orders = get_orders(reffe)
    LagrangianRefFE(typeof(z),p,orders)
  end
  # using map instead of lazy_map seems to kill the compiler.
  ctype_to_reffe = collect(lazy_map(fun,ctype_to_reffe))
  cell_to_reffe = expand_cell_data(ctype_to_reffe,cell_to_ctype)
  cell_to_nodes = get_cell_node_ids(grid)

  cell_to_dofs = _generate_cell_dofs_clagrangian(
    z,
    cell_to_nodes,
    ctype_to_reffe,
    cell_to_ctype,
    glue.node_and_comp_to_dof)

  cell_to_dofs, cell_to_reffe
end

function _generate_cell_dofs_clagrangian(
  z,
  cell_to_nodes,
  ctype_to_reffe,
  cell_to_ctype,
  node_and_comp_to_dof)

  cell_to_dofs = lazy_map(
    Broadcasting(Reindex(node_and_comp_to_dof)),cell_to_nodes)
  Table(cell_to_dofs)
end

function _generate_cell_dofs_clagrangian(
  z::MultiValue,
  cell_to_nodes,
  ctype_to_reffe,
  cell_to_ctype,
  node_and_comp_to_dof)

  ncomps = num_components(z)

  ctype_to_lnode_to_comp_to_ldof = map(get_node_and_comp_to_dof,ctype_to_reffe)
  ctype_to_num_ldofs = map(num_dofs,ctype_to_reffe)

  cell_to_dofs = _allocate_cell_dofs_clagrangian(
    ctype_to_num_ldofs,cell_to_ctype)

  cache = array_cache(cell_to_nodes)
  _fill_cell_dofs_clagrangian!(
    cell_to_dofs,
    cache,
    cell_to_nodes,
    ctype_to_lnode_to_comp_to_ldof,
    cell_to_ctype,
    node_and_comp_to_dof,
    ncomps)

  cell_to_dofs
end

function _allocate_cell_dofs_clagrangian(ctype_to_num_ldofs,cell_to_ctype)
  ncells = length(cell_to_ctype)
  ptrs = zeros(Int32,ncells+1)
  for cell in 1:ncells
    ctype = cell_to_ctype[cell]
    nldofs = ctype_to_num_ldofs[ctype]
    ptrs[cell+1] = nldofs
  end
  length_to_ptrs!(ptrs)
  ndata = ptrs[end]-1
  data = zeros(Int32,ndata)
  Table(data,ptrs)
end

function _fill_cell_dofs_clagrangian!(
  cell_to_dofs,
  cache,
  cell_to_nodes,
  ctype_to_lnode_and_comp_to_ldof,
  cell_to_ctype,
  node_and_comp_to_dof,
  ncomps)

  ncells = length(cell_to_ctype)
  for cell in 1:ncells
    nodes = getindex!(cache,cell_to_nodes,cell)
    ctype = cell_to_ctype[cell]
    lnode_and_comp_to_ldof = ctype_to_lnode_and_comp_to_ldof[ctype]
    p = cell_to_dofs.ptrs[cell]-1
    for (lnode, node) in enumerate(nodes)
      for comp in 1:ncomps
        ldof = lnode_and_comp_to_ldof[lnode][comp]
        dof = node_and_comp_to_dof[node][comp]
        cell_to_dofs.data[p+ldof] = dof
      end
    end
  end
end

