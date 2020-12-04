module FEOperatorsFromWeakForm

using Test
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Algebra
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Integration
using Gridap.Fields
using Gridap.FESpaces
using Gridap.CellData

u(x) = x[1] + x[2]
f(x) = -(3.0*x[1]+x[2]+1.0)
ν(u,x) = (u+1.0)*x[1]
dν(du,x) = x[1]*du

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)

order = 2
reffe = ReferenceFE(:Lagrangian,Float64,order)
V = FESpace(model,reffe,dirichlet_tags="boundary")
U = TrialFESpace(V,u)

degree = 3*order-1
trian = get_triangulation(model)
dΩ = LebesgueMeasure(trian,degree)

a(u,du,v) = ∫( ∇(v)⋅(ν∘(u,identity)*∇(du)) )*dΩ
res(u,v) = a(u,v,u) - ∫(v*f)*dΩ
jac(u,du,v) = a(u,du,v) + ∫( ∇(v)⋅(dν∘(du,identity)*∇(u)) )*dΩ

v = get_cell_shapefuns(V)
du = get_cell_shapefuns_trial(U)
uh = zero(U)

op = FEOperator(res,jac,U,V)

uh = solve(op)
e = u - uh

el2 = sqrt(sum(∫( e*e )*dΩ))
eh1 = sqrt(sum(∫( e*e + ∇(e)⋅∇(e) )*dΩ))

@test el2 < 1.e-8
@test eh1 < 1.e-7

y = ones(num_free_dofs(U))
uh = FEFunction(U,y)
r = residual(op,uh)
A = jacobian(op,uh)
test_fe_operator(op,y,r,≈,jac=A)

#using Gridap.Visualization
#writevtk(trian,"trian",celldata=["e"=>cont])

# Now with autodiff

op = FEOperator(res,U,V)
test_fe_operator(op,y,r,≈,jac=A)

uh = solve(op)
e = u - uh

el2 = sqrt(sum(∫( e*e )*dΩ))
eh1 = sqrt(sum(∫( e*e + ∇(e)⋅∇(e) )*dΩ))

@test el2 < 1.e-8
@test eh1 < 1.e-7

end # module
