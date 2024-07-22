
"""
    struct SimplexifyRefinement <: AdaptivityMethod

Equivalent to `simplexify`, but keeps track of the parent-child relationship between
the original and the refined model.
"""
struct SimplexifyRefinement <: AdaptivityMethod end

function refine(::SimplexifyRefinement,model::UnstructuredDiscreteModel{Dc,Dp};kwargs...) where {Dc,Dp}
  ref_model = simplexify(model;kwargs...)

  polys = get_polytopes(model)
  ctype = get_cell_type(model)
  rrules = expand_cell_data(map(p -> SimplexifyRefinementRule(p;kwargs...),polys),ctype)
  glue = blocked_refinement_glue(rrules)

  return AdaptedDiscreteModel(ref_model,model,glue)
end

struct SimplexifyRefinementRule <: RefinementRuleType end

function SimplexifyRefinementRule(poly::Polytope;kwargs...)
  conn, sp = simplexify(poly;kwargs...)
  coords = get_vertex_coordinates(poly)
  reffes = [LagrangianRefFE(Float64,sp,1)]
  cell_types = fill(1,length(conn))
  ref_grid = UnstructuredDiscreteModel(UnstructuredGrid(coords,Table(conn),reffes,cell_types))
  return RefinementRule(SimplexifyRefinementRule(),poly,ref_grid)
end
