
function FiniteElements(
  domain::DomainStyle,
  model::DiscreteModel,
  basis::Symbol,
  args...;
  kwargs...)

  FiniteElements(domain,model,Val(basis),args...;kwargs...)
end

function FiniteElements(::DomainStyle,model::DiscreteModel,::Val,args...;kwargs...)
  @abstractmethod "The factory function FiniteElements has not been defined for the given arguments"
end

function FiniteElements(
  ::ReferenceDomain,
  model::DiscreteModel,
  basis::Val,
  args...;
  kwargs...)

  cell_reffe = ReferenceFE(model,basis,args...;kwargs...)
  cell_map = get_cell_map(Triangulation(model))
  CellFE(cell_map,cell_reffe)
end

function FiniteElements(
  ::PhysicalDomain,
  model::DiscreteModel,
  basis::Val{:Lagrangian},
  args...;
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
  cell_nodes = lazy_map(evaluate,cell_map,cell_ref_nodes)
  cell_ref_dof_basis = expand_cell_data(ctype_ref_dof_basis,cell_ctype)
  cell_dof_basis = lazy_map(LagrangianDofBasis,cell_ref_dof_basis,cell_nodes)

  # Shape functions
  cell_dof_values = lazy_map(evaluate,cell_dof_basis,cell_prebasis)
  cell_change = lazy_map(inv,cell_dof_values)
  cell_shapefuns = lazy_map(linear_combination,cell_change,cell_prebasis)

  # Build the CellFE
  ctype_num_dofs = map(num_dofs,ctype_reffe)
  ctype_ldof_comp = map(reffe->get_dof_to_comp(reffe),ctype_reffe)
  cell_conformity = CellConformity(cell_reffe)
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

