using Gridap
using Gridap.Geometry, Gridap.ReferenceFEs, Gridap.TensorValues
using Gridap.Fields, Gridap.FESpaces, Gridap.CellData

using FillArrays

using Meshes: viz
using CairoMakie

# TODO: Add support for GridapMakie, which is what Meshes.jl uses for plotting underneath
# using GridapMakie
# using Makie
# using GridapMakie: PlotGridMesh, PlotGrid
# Makie.plottype(::Grid) = PlotGridMesh
# function Makie.convert_arguments(::Type{PlotGridMesh}, grid::Grid)
#   (PlotGrid(grid), )
# end
# function Makie.convert_arguments(t::Type{<:Union{Makie.Wireframe, Makie.Scatter}}, grid::Grid)
#   Makie.convert_arguments(t, PlotGrid(grid))
# end
# Makie.args_preferred_axis(pg::PlotGrid)= num_point_dims(pg.grid)<=2 ? Makie.Axis : Makie.LScene

model = CartesianDiscreteModel((0,1,0,1),(2,2))

pmodel = Gridap.Geometry.PolytopalDiscreteModel(model)
vmodel = Gridap.Geometry.voronoi(Gridap.Geometry.simplexify(model))
polys = get_polytopes(vmodel)

Geometry.restrict(vmodel,[1,2,3,4])

writevtk(vmodel,"tmp/polygonal_model")
viz(vmodel;color=1:num_cells(vmodel),showpoints=true,showsegments=true)

Ω = Triangulation(vmodel)
Γ = Boundary(vmodel)
Λ = Skeleton(vmodel)

order = 1
u_exact(x) = x[1]^order + x[2]^order

dΩ = Measure(Ω,2*order)
dΓ = Measure(Γ,2*order)
dΛ = Measure(Λ,2*order)

trian = get_triangulation(Ω)
grid  = get_grid(trian) 

sum(∫(1)dΩ)

V = FESpaces.PolytopalFESpace(Ω,Float64,order,space=:P)

# L2 projection
mass_lhs(u,v) = ∫(u⋅v)dΩ
mass_rhs(v) = ∫(v*u_exact)dΩ

op = AffineFEOperator(mass_lhs,mass_rhs,V,V)
uh = solve(op)

eh = uh - u_exact
sum(∫(eh⋅eh)dΩ)

writevtk(
  Triangulation(vmodel),"tmp/polygonal_trian",
  cellfields = Dict("uh" => uh, "u_exact" => u_exact, "eh" => eh),
)

# DG Poisson
β = 100.0
f(x) = -Δ(u_exact)(x)
nΛ = get_normal_vector(Λ)
βΛ = CellField(β ./ get_array(∫(1)dΛ),Λ)
nΓ = get_normal_vector(Γ)
βΓ = CellField(β ./ get_array(∫(1)dΓ),Γ)

lap_lhs(u,v) = ∫(∇(u)⋅∇(v))dΩ + ∫(βΛ*jump(u*nΛ)⋅jump(v*nΛ) - mean(∇(u))⋅jump(v*nΛ) - mean(∇(v))⋅jump(u*nΛ))dΛ + ∫(βΓ*u*v - (∇(u)⋅nΓ)*v - (∇(v)⋅nΓ)*u)dΓ
lap_rhs(v) = ∫(v*f)dΩ + ∫(βΓ*u_exact*v - (∇(v)⋅nΓ)*u_exact)dΓ

op = AffineFEOperator(lap_lhs,lap_rhs,V,V)
uh = solve(op)

eh = uh - u_exact
sum(∫(eh⋅eh)dΩ)

###########################
# 3D

model = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))

pmodel = Gridap.Geometry.PolytopalDiscreteModel(model)
polys = get_polytopes(pmodel)
writevtk(pmodel,"tmp/polygonal_model";append=false)
