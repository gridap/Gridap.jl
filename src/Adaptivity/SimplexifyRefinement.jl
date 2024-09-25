
"""
    struct SimplexifyRefinement <: AdaptivityMethod

Equivalent to `simplexify`, but keeps track of the parent-child relationship between
the original and the refined model.
"""
struct SimplexifyRefinement <: AdaptivityMethod end

function refine(::SimplexifyRefinement,model::CartesianDiscreteModel{Dc,Dp}) where {Dc,Dp}
  _simplex_refine(model)
end

function refine(::SimplexifyRefinement,model::UnstructuredDiscreteModel{Dc,Dp}) where {Dc,Dp}
  _simplex_refine(model)
end

function _simplex_refine(model::DiscreteModel{Dc,Dp}) where {Dc,Dp}
  ref_model = simplexify(model)

  polys = get_polytopes(model)
  ctype = get_cell_type(model)
  rrules = expand_cell_data(map(p -> SimplexifyRefinementRule(p),polys),ctype)
  glue = blocked_refinement_glue(rrules)

  return AdaptedDiscreteModel(ref_model,model,glue)
end

struct SimplexifyRefinementRule <: RefinementRuleType end

function SimplexifyRefinementRule(poly::Polytope)
  conn, sp = simplexify(poly)
  coords = get_vertex_coordinates(poly)
  reffes = [LagrangianRefFE(Float64,sp,1)]
  cell_types = fill(1,length(conn))
  ref_grid = UnstructuredDiscreteModel(UnstructuredGrid(coords,Table(conn),reffes,cell_types))
  return RefinementRule(SimplexifyRefinementRule(),poly,ref_grid)
end
