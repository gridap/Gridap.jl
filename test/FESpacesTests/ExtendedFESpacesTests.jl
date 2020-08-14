module ExtendedFESpacesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Algebra
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.FESpaces: ExtendedVector
using Gridap.Integration
using Gridap.Fields


oldcell_to_cell = [1,2,-1,3,-2,-3,-4]
cell_to_oldcell = [1,2,4]
void_to_oldcell = [3,5,6,7]
void_to_val = Float64[-10,-20,-30,-40]
cell_to_val = Float64[10,20,30]
a = ExtendedVector(
  void_to_val,
  cell_to_val,
  oldcell_to_cell,
  void_to_oldcell,
  cell_to_oldcell)

r = [10.0, 20.0, -10.0, 30.0, -20.0, -30.0, -40.0]
test_array(a,r)


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
oldcell_to_is_in = collect1d(apply(is_in,oldcell_to_coods))

incell_to_cell = findall(oldcell_to_is_in)
outcell_to_cell = findall(collect(Bool, .! oldcell_to_is_in))

model_in = DiscreteModel(model,incell_to_cell)

V = TestFESpace(
  model=model_in,valuetype=Float64,reffe=:Lagrangian,
  order=1,conformity=:H1,dof_space=:physical)
@test isa(V,ExtendedFESpace)



model_out = DiscreteModel(model,outcell_to_cell)

trian_in = RestrictedTriangulation(trian, incell_to_cell)

trian_Γ = InterfaceTriangulation(model_in,model_out)

quad_in = CellQuadrature(trian_in,2*order)

quad_Γ = CellQuadrature(trian_Γ,2*order)

reffes = [LagrangianRefFE(VectorValue{2,Float64},get_polytope(p),order) for p in get_reffes(trian_in)]

V_in = DiscontinuousFESpace(reffes,trian_in)

V = ExtendedFESpace(V_in, trian_in)

test_single_field_fe_space(V)

U = TrialFESpace(V)
test_single_field_fe_space(U)

u(x) = VectorValue(x[1]+x[2], x[1])

uh = interpolate(u,U)


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

end # module
