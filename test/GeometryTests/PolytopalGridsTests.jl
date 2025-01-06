using Gridap
using Gridap.Geometry, Gridap.ReferenceFEs, Gridap.TensorValues
using Gridap.Fields, Gridap.FESpaces, Gridap.CellData

using FillArrays

using Meshes: viz
using CairoMakie

model = CartesianDiscreteModel((0,1,0,1),(4,4))

pmodel = Geometry.PolytopalDiscreteModel(model)
vmodel = Geometry.voronoi(simplexify(model))
polys = get_polytopes(vmodel)
# viz(vmodel;color=1:num_cells(vmodel),showpoints=true,showsegments=true)

Ω = Triangulation(vmodel)
dΩ = Measure(Ω,2)

V = FESpaces.PolytopalFESpace(Ω,1)

u_exact(x) = x[1] + x[2]
a(u,v) = ∫(u⋅v)dΩ
l(v) = ∫(v*u_exact)dΩ

op = AffineFEOperator(a,l,V,V)
uh = solve(op)

eh = uh - u_exact
sum(∫(eh⋅eh)dΩ)
