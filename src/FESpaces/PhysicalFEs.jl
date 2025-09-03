
"""
"""
function FiniteElements(::DomainStyle,model::DiscreteModel,basis::ReferenceFEName,args...;kwargs...)
  @abstractmethod "The factory function FiniteElements has not been defined for the given arguments"
end

function FiniteElements(
  ::ReferenceDomain,
  model::DiscreteModel,
  basis::ReferenceFEName,
  args...;
  conformity=nothing,
  kwargs...)

  cell_reffe = ReferenceFE(model,basis,args...;kwargs...)
  conf = Conformity(testitem(cell_reffe),conformity)
  CellFE(model,cell_reffe,conf)
end

function FiniteElements(
  ::PhysicalDomain,
  model::DiscreteModel,
  basis::Lagrangian,
  args...;
  conformity=nothing,
  kwargs...)

  # Reference FEs and cell_map
  cell_reffe = ReferenceFE(model,basis,args...;kwargs...)
  cell_map = get_cell_map(Triangulation(model))
  ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)

  # prebasis # TODO we need a better conditioned cell_prebasis
  ctype_prebasis = map(get_prebasis,ctype_reffe)
  cell_prebasis = expand_cell_data(ctype_prebasis,cell_ctype)

  # dof basis
  ctype_ref_dof_basis = map(get_dof_basis,ctype_reffe)
  ctype_ref_nodes = map(i->i.nodes,ctype_ref_dof_basis)
  cell_ref_nodes = expand_cell_data(ctype_ref_nodes,cell_ctype)
  cell_node_ids = lazy_map(evaluate,cell_map,cell_ref_nodes)
  cell_ref_dof_basis = expand_cell_data(ctype_ref_dof_basis,cell_ctype)
  cell_dof_basis = lazy_map(LagrangianDofBasis,cell_ref_dof_basis,cell_node_ids)

  # Shape functions
  cell_dof_values = lazy_map(evaluate,cell_dof_basis,cell_prebasis)
  cell_change = lazy_map(inv,cell_dof_values)
  cell_shapefuns = lazy_map(linear_combination,cell_change,cell_prebasis)

  # Build the CellFE
  ctype_num_dofs = map(num_dofs,ctype_reffe)
  ctype_ldof_comp = map(reffe->get_dof_to_comp(reffe),ctype_reffe)
  conf = Conformity(testitem(cell_reffe),conformity)
  cell_conformity = CellConformity(cell_reffe,conf)
  cell_shapefuns_domain = PhysicalDomain()
  cell_dof_basis_domain = cell_shapefuns_domain
  max_order = maximum(map(get_order,ctype_reffe))

  CellFE(
    cell_ctype,
    ctype_num_dofs,
    ctype_ldof_comp,
    cell_conformity,
    cell_shapefuns,
    cell_dof_basis,
    cell_shapefuns_domain,
    cell_dof_basis_domain,
    max_order)
end
