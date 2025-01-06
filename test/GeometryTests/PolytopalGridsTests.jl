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
Γ = Boundary(vmodel)
Λ = Skeleton(vmodel)

order = 2
u_exact(x) = x[1]^order + x[2]^order

dΩ = Measure(Ω,2*order)
dΓ = Measure(Γ,2*order)

sum(∫(1)dΩ)

V = FESpaces.PolytopalFESpace(Ω,order)

# L2 projection

mass_lhs(u,v) = ∫(u⋅v)dΩ
mass_rhs(v) = ∫(v*u_exact)dΩ

op = AffineFEOperator(mass_lhs,mass_rhs,V,V)
uh = solve(op)

eh = uh - u_exact
sum(∫(eh⋅eh)dΩ)

# DG Poisson
f(x) = Δ(u_exact)(x)
lap_lhs(u,v) = ∫(∇(u)⋅∇(v))dΩ
lap_rhs(v) = ∫(v*f)dΩ
