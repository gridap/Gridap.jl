using Gridap
using Gridap.Geometry

using Meshes, CairoMakie

model = CartesianDiscreteModel((0,1,0,1),(4,4))

pmodel = Geometry.PolytopalDiscreteModel(model)
vmodel = Geometry.voronoi(simplexify(model))
viz(vmodel;color=1:num_cells(vmodel),showpoints=true,showsegments=true)

polys = get_polytopes(vmodel)
