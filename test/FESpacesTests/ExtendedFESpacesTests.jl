module ExtendedFESpacesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Algebra
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Integration
using Gridap.Fields
using Gridap.CellData

n = 10
mesh = (n,n)
domain = (0,1,0,1) .- 1
order = 1
model = CartesianDiscreteModel(domain, mesh)

trian = Triangulation(model)

const R = 0.7

function is_in(coords)
  n = length(coords)
  x = (1/n)*sum(coords)
  d = x[1]^2 + x[2]^2 - R^2
  d < 0
end

oldcell_to_coods = get_cell_coordinates(trian)
oldcell_to_is_in = collect1d(lazy_map(is_in,oldcell_to_coods))

model_in = DiscreteModel(model,oldcell_to_is_in)
@test isa(model_in,RestrictedDiscreteModel)

@test trian === get_triangulation(get_parent_model(model_in))

V = TestFESpace(model_in,ReferenceFE(:Lagrangian,Float64,order),conformity=:H1)
@test isa(V,FESpaces.ExtendedFESpace)
test_single_field_fe_space(V)

U = TrialFESpace(V)
test_single_field_fe_space(U)

oldcell_to_is_out = lazy_map(!,oldcell_to_is_in)
model_out = DiscreteModel(model,oldcell_to_is_out)

Ω = Triangulation(model)
Ω_in = Triangulation(model_in)
Ω_out = Triangulation(model_out)
Γ = InterfaceTriangulation(model_in,model_out)

degree = 2
dΩ = LebesgueMeasure(Ω,degree)
dΩ_in = LebesgueMeasure(Ω_in,degree)
dΩ_out = LebesgueMeasure(Ω_out,degree)
dΓ = LebesgueMeasure(Γ,degree)

u(x) = x[1]+1

dv = get_cell_shapefuns(V)
du = get_cell_shapefuns_trial(U)
uh = interpolate(u,U)

x = get_cell_points(Ω)
x_in = get_cell_points(Ω_in)
x_out = get_cell_points(Ω_out)
x_Γ = get_cell_points(Γ)

r = dv(x)
test_array(r,collect(r))
r = dv(x_in)
test_array(r,collect(r))
r = dv(x_out)
test_array(r,collect(r))
r = dv.⁺(x_Γ)
test_array(r,collect(r))

r = ∇(dv)(x)
test_array(r,collect(r))
r = ∇(dv)(x_in)
test_array(r,collect(r))
r = ∇(dv)(x_out)
test_array(r,collect(r))
r = ∇(dv).⁺(x_Γ)
test_array(r,collect(r))

r = du(x)
test_array(r,collect(r))
r = du(x_in)
test_array(r,collect(r))
r = du(x_out)
test_array(r,collect(r))
r = du.⁺(x_Γ)
test_array(r,collect(r))

r = ∇(du)(x)
test_array(r,collect(r))
r = ∇(du)(x_in)
test_array(r,collect(r))
r = ∇(du)(x_out)
test_array(r,collect(r))
r = ∇(du).⁺(x_Γ)
test_array(r,collect(r))

r = uh(x)
test_array(r,collect(r))
r = uh(x_in)
test_array(r,collect(r))
r = uh(x_out)
test_array(r,collect(r))
r = uh.⁺(x_Γ)
test_array(r,collect(r))

r = ∇(uh)(x)
test_array(r,collect(r))
r = ∇(uh)(x_in)
test_array(r,collect(r))
r = ∇(uh)(x_out)
test_array(r,collect(r))
r = ∇(uh).⁺(x_Γ)
test_array(r,collect(r))

#using Gridap.Visualization
#writevtk(Ω,"Omega",cellfields=["uh"=>uh])
#writevtk(Ω_in,"Omega_in",cellfields=["uh"=>uh])
#writevtk(Ω_out,"Omega_out",cellfields=["uh"=>uh])
#writevtk(Γ,"gamma",cellfields=["uh"=>uh.plus])

a(u,v) =
  ∫(v⊙u)*dΩ #+
  #∫(v⊙u)*dΩ_in +
  #∫(v⊙u)*dΩ_out# +
  #∫(jump(v)⊙jump(u) + jump(∇(v))⊙jump(∇(u)))*dΓ

l(v) =
  ∫(v⊙u)*dΩ #+
  #∫(v⊙u)*dΩ_in +
  #∫(v⊙u)*dΩ_out# +
  #∫(jump(v))*dΓ

op = AffineFEOperator(a,l,U,V)

kk




uh_in = restrict(uh,trian_in)

uh_Γ = restrict(uh,trian_Γ)

t_in = AffineFETerm( (u,v) -> v⊙u, (v) -> v⊙u, trian_in, quad_in)
op_in = AffineFEOperator(U,V,t_in)

quad = CellQuadrature(trian,2*order)

t_Ω = AffineFETerm( (u,v) -> v⊙u, (v) -> v⊙u, trian, quad)
op_Ω = AffineFEOperator(U,V,t_Ω)

@test get_vector(op_in) ≈ get_vector(op_Ω)

t_Γ = AffineFETerm( (u,v) -> jump(v)⊙jump(u) + jump(ε(v))⊙jump(ε(u)), (v) -> jump(v)⊙u, trian_Γ, quad_Γ)
op_Γ = AffineFEOperator(U,V,t_Γ)

q_in = get_coordinates(quad_in)
collect(evaluate(uh_in,q_in))

q = get_coordinates(quad)
collect(evaluate(uh,q))

q_Γ = get_coordinates(quad_Γ)
collect(evaluate(jump(uh_Γ),q_Γ))
collect(evaluate(jump(∇(uh_Γ)),q_Γ))
collect(evaluate(jump(ε(uh_Γ)),q_Γ))

V = TestFESpace(model=model_in,valuetype=Float64,reffe=:Lagrangian,order=2,conformity=:H1)
@test isa(V,ExtendedFESpace)

V = TestFESpace(model=model_in,valuetype=Float64,reffe=:Lagrangian,order=2,conformity=:H1,constraint=:zeromean)
uh = FEFunction(V,rand(num_free_dofs(V)))
uh_in = restrict(uh,trian_in)
@test sum(integrate(uh_in,trian_in,quad_in)) + 1 ≈ 1

V = TestFESpace(model=model,valuetype=Float64,reffe=:Lagrangian,order=2,conformity=:H1)
@test !isa(V,ExtendedFESpace)

V = TestFESpace(triangulation=trian_in,valuetype=Float64,reffe=:Lagrangian,order=2,conformity=:L2)
@test isa(V,ExtendedFESpace)

V = TestFESpace(triangulation=trian_in,valuetype=Float64,reffe=:Lagrangian,order=2,conformity=:L2,constraint=:zeromean)
uh = FEFunction(V,rand(num_free_dofs(V)))
uh_in = restrict(uh,trian_in)
@test sum(integrate(uh_in,trian_in,quad_in)) + 1 ≈ 1

V = TestFESpace(triangulation=trian,valuetype=Float64,reffe=:Lagrangian,order=2,conformity=:L2)
@test !isa(V,ExtendedFESpace)

#using Gridap.Visualization
#writevtk(trian,"trian",cellfields=["uh"=>uh])
#writevtk(trian_in,"trian_in",cellfields=["uh"=>uh_in])

V_in = TestFESpace(model=model_in,valuetype=Float64,reffe=:Lagrangian,order=2,conformity=:H1)
V = TestFESpace(model=model,valuetype=Float64,reffe=:Lagrangian,order=2,conformity=:H1)

vh_in = interpolate(V_in) do x
    x[1]
end
vh_in = interpolate(vh_in, V_in)
vh = interpolate(vh_in, V)

#using Gridap.Visualization
#writevtk(trian,"trian",cellfields=["vh"=>vh,"vh_in"=>vh_in])


V_in = TestFESpace(
  model=model_in,valuetype=Float64,reffe=:Lagrangian,
  order=1,conformity=:H1,dof_space=:physical)

V_out = TestFESpace(
  model=model_out,valuetype=Float64,reffe=:Lagrangian,
  order=2,conformity=:H1,dof_space=:physical)

cell_axes = Fields.create_array_of_blocked_axes(get_cell_axes(V_in),get_cell_axes(V_out))
#using BenchmarkTools
#cache = array_cache(cell_axes)
#@btime getindex!($cache,$cell_axes,1)

end # module
