
# The user is not expected to call this function. Use the factory function instead
function _ConformingFESpace(
  model::DiscreteModel,
  face_labeling::FaceLabeling,
  cell_reffe::AbstractVector{<:ReferenceFE},
  domain_style::DomainStyle,
  conformity::Conformity,
  dirichlet_tags,
  dirichlet_components)

  if domain_style == ReferenceDomain()
    cell_shapefuns, cell_dof_basis = _compute_reference_cell_space(cell_reffe,cell_map)
  else
    cell_shapefuns, cell_dof_basis = _compute_physical_cell_space(cell_reffe,cell_map)
  end

end

function _compute_reference_cell_space(cell_reffe,cell_map)
  # TODO only for Lagrangian at this moment
  ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
  ctype_shapefuns = map(get_shapefuns,ctype_reffe)
  ctype_dof_basis = map(get_dof_basis,ctype_reffe)
  cell_shapefuns = expand_cell_data(ctype_shapefuns,cell_ctype)
  cell_dof_basis = expand_cell_data(ctype_dof_basis,cell_ctype)
  cell_shapefuns, cell_dof_basis
end

function _compute_physical_cell_space(cell_reffe,cell_map)
  # For any reffe whose dof basis implements the DofBasisMap
  ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
  ctype_prebasis = map(get_prebasis,ctype_reffe)
  ctype_ref_dof_basis = map(get_dof_basis,ctype_reffe)
  cell_prebasis = expand_cell_data(ctype_prebasis,cell_ctype)
  cell_ref_dof_basis = expand_cell_data(ctype_ref_dof_basis,cell_ctype)
  cell_dof_basis = lazy_map(DofBasisMap(),cell_ref_dof_basis,cell_map)
  cell_dof_values = lazy_map(evaluate,cell_dof_basis,cell_prebasis)
  cell_change = lazy_map(inv,cell_dof_values)
  cell_shapefuns = lazy_map(linear_combination,cell_change,cell_prebasis)
  cell_shapefuns, cell_dof_basis
end
