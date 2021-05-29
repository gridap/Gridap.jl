module MultiFieldFEOperatorsTests

using Gridap.Arrays
using Gridap.Algebra
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
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
dΩ = Measure(Ω,degree)

sdegree = order
Γ = SkeletonTriangulation(model)
dΓ = Measure(Γ,sdegree)

V = TestFESpace(model,ReferenceFE(lagrangian,Float64,order);conformity=:H1)
Q = TestFESpace(model,ReferenceFE(lagrangian,Float64,order-1),conformity=:L2)

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
test_fe_operator(op,get_free_dof_values(xh),b)

r((u,p),(v,q)) = ∫( v*(u*u) + v*p*u - q*p - v*4 + q )*dΩ
j((u,p),(du,dp),(v,q)) = ∫(2*v*u*du + v*dp*u + v*p*du - q*dp)*dΩ

op = FEOperator(r,j,X,Y)
xh = zero(X)
b = residual(op,xh)
A = jacobian(op,xh)
test_fe_operator(op,get_free_dof_values(xh),b)

@test_broken begin
op_auto = FEOperator(r,X,Y)
A_auto = jacobian(op_auto,xh)
@test A ≈ A_auto
true
end

r_const((u,p),(v,q)) = -1.0 * (∫( v*1.0 )*dΩ)
op_const = FEOperator(r_const,j,X,Y)
b_const = residual(op_const,xh)
test_fe_operator(op_const,get_free_dof_values(xh),b_const)

end # module
