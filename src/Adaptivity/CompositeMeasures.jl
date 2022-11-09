
"""
  Composite Measure

  Measure such that the integration and target triangulations are different. 

  - ttrian: Target triangulation, where the domain contribution lives.
  - itrian: Integration triangulation, where the integration takes place.
  - quad  : CellQuadrature, defined in itrian
"""
struct CompositeMeasure{A<:Triangulation,B<:Triangulation,C<:CellQuadrature} <: GridapType
  ttrian :: A
  itrian :: B
  quad   :: C
end

CellData.Measure(ttrian::Triangulation,itrian::Triangulation,args...;kwargs...) = CompositeMeasure(ttrian,itrian,args...;kwargs...)

function CompositeMeasure(ttrian::Triangulation,itrian::Triangulation,args...;kwargs...)
  quad = CellQuadrature(itrian,args...;kwargs...)
  return CompositeMeasure(ttrian,itrian,quad)
end

get_cell_points(a::CompositeMeasure) = get_cell_points(a.quad)

function CellData.integrate(f,b::CompositeMeasure)
  ic = CellData.integrate(f,b.quad)
  cont = DomainContribution()
  tc = move_contributions(ic,b.itrian,b.ttrian)
  add_contribution!(cont,b.ttrian,tc)
  cont
end

function (*)(a::CellData.Integrand,b::CompositeMeasure)
  CellData.integrate(a.object,b)
end

(*)(b::CompositeMeasure,a::CellData.Integrand) = a*b