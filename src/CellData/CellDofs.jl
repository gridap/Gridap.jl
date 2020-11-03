
"""
"""
struct CellDof{DS} <: CellDatum
  cell_dof::AbstractArray{<:Union{Dof,AbstractArray{<:Dof}}}
  trian::Triangulation
  domain_style::DS
end

get_cell_data(f::CellDof) = f.cell_dof
get_triangulation(f::CellDof) = f.trian
DomainStyle(::Type{CellDof{DS}}) where DS = DS()

function change_domain(a::CellDof,::ReferenceDomain,::PhysicalDomain)
  @notimplemented
end

function change_domain(a::CellDof,::PhysicalDomain,::ReferenceDomain)
  @notimplemented
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
  lazy_map(evaluate,get_cell_data(s),get_cell_data(b))
end
