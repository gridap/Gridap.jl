module MultiFieldFEOperatorsTests

using Gridap.Arrays
using Gridap.Algebra
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
using Gridap.Integration
using Gridap.MultiField
using Gridap.CellData
using SparseArrays
using Test
using LinearAlgebra
using Gridap.ReferenceFEs

using Gridap.MultiField

order = 2

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

degree = order
Ω = get_triangulation(model)
dΩ = LebesgueMeasure(Ω,degree)

sdegree = order
Γ = SkeletonTriangulation(model)
dΓ = LebesgueMeasure(Γ,sdegree)

V = TestFESpace(model,ReferenceFE(:Lagrangian,Float64,order);conformity=:H1)
Q = TestFESpace(model,ReferenceFE(:Lagrangian,Float64,order-1),conformity=:L2)

U = TrialFESpace(V)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

a((u,p),(v,q)) =
  ∫( v*u + v*p - q*p )*dΩ +
  ∫( jump(v)*mean(u) + jump(∇(q))⋅jump(∇(p)) )*dΓ

l((v,q)) = ∫( v*4 + q )*dΩ

op = AffineFEOperator(a,l,X,Y)

xh = zero(X)
b = residual(op,xh)
A = jacobian(op,xh)
test_fe_operator(op,get_free_values(xh),b)


r((u,p),(v,q)) = ∫( v*(u*u) + v*p*u - q*p - v*4 + q )*dΩ

function j(x,dx,y)
  u,p = x
  du,dp = dx
  v,q = y
  2*v*u*du + v*dp*u + v*p*du - q*dp
end


op = FEOperator(X,Y,t_Ω_auto)
xh = zero(X)
b = residual(op,xh)
A = jacobian(op,xh)
test_fe_operator(op,get_free_values(xh),b)

t_Ω = FETerm(r,j,trian,quad)
t_Ω_auto = FETerm(r,trian,quad)

x = FEFunction(X,rand(num_free_dofs(X)))
dx = get_cell_basis(X)
y = get_cell_basis(Y)

cell_r = get_cell_residual(t_Ω,x,y)
cell_j = get_cell_jacobian(t_Ω,x,dx,y)

cell_r_auto = get_cell_residual(t_Ω_auto,x,y)
cell_j_auto = get_cell_jacobian(t_Ω_auto,x,dx,y)
test_array(cell_r_auto,cell_r,≈)
test_array(cell_j_auto,cell_j,≈)



end # module
