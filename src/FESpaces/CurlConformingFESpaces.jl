
# @santiagobadia : Probably not needed, the same as div...
"""
    CurlConformingFESpace(
      reffes::Vector{<:ReferenceFE},
      model::DiscreteModel,
      face_labeling::FaceLabeling,
      dirichlet_tags)
"""
function CurlConformingFESpace(
  reffes::Vector{<:ReferenceFE},
  model::DiscreteModel,
  face_labeling::FaceLabeling,
  dirichlet_tags,
  is_ref)

  grid_topology = get_grid_topology(model)

  cell_dofs, nfree, ndirichlet, dirichlet_dof_tag, dirichlet_cells = compute_conforming_cell_dofs(
    reffes,grid_topology,face_labeling,dirichlet_tags)

  ntags = length(dirichlet_tags)

  grid = get_grid(model)
  cell_to_ctype = get_cell_type(grid_topology)
  cell_map = get_cell_map(grid)

  cell_shapefuns, cell_dof_basis = compute_cell_space(reffes,cell_to_ctype,cell_map,Val(is_ref))

  UnconstrainedFESpace(
    nfree,
    ndirichlet,
    cell_dofs,
    cell_shapefuns,
    cell_dof_basis,
    cell_map,
    dirichlet_dof_tag,
    dirichlet_cells,
    ntags)

end

function _compute_hcurl_cell_space(reffes, cell_to_ctype, cell_map)
  #TODO: fine for structured hex meshes, but not otherwise
  # compute_cell_space(reffes,cell_to_ctype,cell_map)
end
