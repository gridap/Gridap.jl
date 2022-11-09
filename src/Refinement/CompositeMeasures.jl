
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

function Geometry.move_contributions(scell_to_val::AbstractArray, strian::RefinedTriangulation, ttrian::Triangulation)
  smodel = get_refined_model(strian)
  @check get_parent(smodel) === get_background_model(ttrian)

  tcell_to_val = move_contributions(scell_to_val,get_refinement_glue(smodel))
  return tcell_to_val
end

function Geometry.move_contributions(scell_to_val::AbstractArray, glue::RefinementGlue)
  tcell_to_scells = glue.c2f_faces_map
  k = Geometry.CombineContributionsMap(scell_to_val)
  tcell_to_val = lazy_map(k,tcell_to_scells)
  return tcell_to_val
end