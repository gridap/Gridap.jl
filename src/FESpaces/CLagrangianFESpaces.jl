
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
    CLagrangianFESpace(::Type{T},grid::Triangulation) where T
"""
function CLagrangianFESpace(::Type{T},grid::Triangulation,trian::Triangulation=grid) where T
  vector_type = Vector{_dof_type(T)}
  node_to_tag = fill(Int8(UNSET),num_nodes(grid))
  tag_to_mask = fill(_default_mask(T),0)
  CLagrangianFESpace(T,grid,vector_type,node_to_tag,tag_to_mask,trian)
end

"""
    CLagrangianFESpace(
    ::Type{T},
    grid::Triangulation,
    vector_type::Type,
    node_to_tag::AbstractVector,
    tag_to_mask::AbstractVector) where T
"""
function CLagrangianFESpace(
  ::Type{T},
  grid::Triangulation,
  vector_type::Type,
  node_to_tag::AbstractVector{<:Integer},
  tag_to_masks::AbstractVector,
  trian::Triangulation=grid
  ) where T

  z = zero(T)
  glue, dirichlet_dof_tag = _generate_node_to_dof_glue_component_major(
    z,node_to_tag,tag_to_masks)
  cell_reffe = _generate_cell_reffe_clagrangian(z,grid)
  cell_dofs_ids = _generate_cell_dofs_clagrangian(z,grid,glue,cell_reffe)
  cell_shapefuns = lazy_map(get_shapefuns,cell_reffe)
  cell_dof_basis = lazy_map(get_dof_basis,cell_reffe)
  rd = ReferenceDomain()
  fe_basis, fe_dof_basis = compute_cell_space(
    cell_shapefuns,cell_dof_basis,rd,rd,trian)
  cell_is_dirichlet = _generate_cell_is_dirichlet(cell_dofs_ids)

  nfree = length(glue.free_dof_to_node)
  ndirichlet = length(glue.dirichlet_dof_to_node)
  ntags = length(tag_to_masks)
  dirichlet_cells = collect(Int32,findall(cell_is_dirichlet))

  space = UnconstrainedFESpace(
    vector_type,
    nfree,
    ndirichlet,
    cell_dofs_ids,
    fe_basis,
    fe_dof_basis,
    cell_is_dirichlet,
    dirichlet_dof_tag,
    dirichlet_cells,
    ntags,
    glue)

  space
end

function _use_clagrangian(trian,cell_reffe,conf)
  false
end

function _use_clagrangian(trian::Triangulation,cell_reffe,conf::H1Conformity)
  ctype_reffe1, cell_ctype1 = compress_cell_data(cell_reffe)
  ctype_reffe2 = get_reffes(trian)
  if length(ctype_reffe1) != 1 || length(ctype_reffe2) != 1
    return false
  end
  reffe1 = first(ctype_reffe1)
  reffe2 = first(ctype_reffe2)
  if get_orders(reffe1) != get_orders(reffe2)
    return false
  end
  if get_order(reffe1) != 1 # This can be relaxed in the future
    return false
  end
  return true
end

function _unsafe_clagrangian(
  cell_reffe,
  grid,
  labels,
  vector_type,
  dirichlet_tags,
  dirichlet_masks,
  trian=grid)

  ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
  prebasis = get_prebasis(first(ctype_reffe))
  T = return_type(prebasis)
  # Next line assumes linear grids
  node_to_tag = get_face_tag_index(labels,dirichlet_tags,0)
  _vector_type = vector_type === nothing ? Vector{Float64} : vector_type
  tag_to_mask = dirichlet_masks === nothing ? fill(_default_mask(T),length(dirichlet_tags)) : dirichlet_masks
  CLagrangianFESpace(T,grid,_vector_type,node_to_tag,tag_to_mask,trian)
end

# Helpers

_default_mask(::Type) = true
_default_mask(::Type{T}) where T <: MultiValue = ntuple(i->true,Val{length(T)}())

_dof_type(::Type{T}) where T = T
_dof_type(::Type{T}) where T<:MultiValue = eltype(T)

function _generate_cell_is_dirichlet(cell_dofs)
  cell_is_dirichlet = Vector{Bool}(undef,length(cell_dofs))
  cache = array_cache(cell_dofs)
  isnegative(dof) = dof<0
  for cell in 1:length(cell_dofs)
    dofs = getindex!(cache,cell_dofs,cell)
    cell_is_dirichlet[cell] = any(isnegative,dofs)
  end
  cell_is_dirichlet
end

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
  @check length(testitem(tag_to_masks)) == ncomps
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

function _generate_cell_reffe_clagrangian(z,grid)
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
end

function _generate_cell_dofs_clagrangian(z,grid,glue,cell_to_reffe)
  ctype_to_reffe, cell_to_ctype = compress_cell_data(cell_to_reffe)
  cell_to_nodes = get_cell_node_ids(grid)
  cell_to_dofs = _generate_cell_dofs_clagrangian(
    z,
    cell_to_nodes,
    ctype_to_reffe,
    cell_to_ctype,
    glue.node_and_comp_to_dof)
  cell_to_dofs
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

