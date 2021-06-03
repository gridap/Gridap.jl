
"""
"""
struct CellDof{DS} <: CellDatum
  cell_dof::AbstractArray{<:Union{Dof,AbstractArray{<:Dof}}}
  trian::Triangulation
  domain_style::DS
end

get_data(f::CellDof) = f.cell_dof
get_triangulation(f::CellDof) = f.trian
DomainStyle(::Type{CellDof{DS}}) where DS = DS()

function change_domain(a::CellDof,::ReferenceDomain,::PhysicalDomain)
  @notimplemented
end

function change_domain(a::CellDof,::PhysicalDomain,::ReferenceDomain)
  @notimplemented
end

function change_domain(a::CellDof,target_trian::Triangulation,target_domain::DomainStyle)
  @notimplemented
end

function change_domain(a::CellDof,target_trian::RestrictedTriangulation,target_domain::DomainStyle)
  @notimplementedif DomainStyle(a) != target_domain
  trian_a = get_triangulation(a)
  @notimplementedif ! have_compatible_domains(trian_a,get_background_triangulation(target_trian))
  cell_dof = get_data(a)
  tcell_to_cell = get_cell_to_bgcell(target_trian)
  tcell_dof = lazy_map(Reindex(cell_dof),tcell_to_cell)
  CellDof(tcell_dof,target_trian,DomainStyle(a))
end

function get_cell_points(dofs::CellDof)
  cell_to_x = lazy_map(get_nodes,get_data(dofs))
  CellPoint(cell_to_x,get_triangulation(dofs),DomainStyle(dofs))
end

# Evaluation of CellDof

(a::CellDof)(f) = evaluate(a,f)

function evaluate!(cache,s::CellDof,f::CellField)

  trian_f = get_triangulation(f)
  trian_s = get_triangulation(s)

  if trian_f !== trian_s
    @unreachable """\n
    A CellDof can only be evaluated on a CellField defined on the same Triangulation.
    """
  end

  b = change_domain(f,s.domain_style)
  lazy_map(evaluate,get_data(s),get_data(b))
end

function evaluate!(cache, ::CellField, ::CellDof)
  @unreachable """\n
  CellField (f) objects cannot be evaluated at CellDof (s) objects.
  However, CellDofs objects can be evaluated at CellField objects.
  Did you mean evaluate(f,s) instead of evaluate(s,f), i.e.
  f(s) instead of s(f)?
  """
end
