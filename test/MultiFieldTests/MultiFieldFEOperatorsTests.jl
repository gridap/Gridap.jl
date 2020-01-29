module MultiFieldFEOperatorsTests

using Gridap.Arrays
using Gridap.Algebra
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
using Gridap.Integration
using SparseArrays
using Test

using Gridap.MultiField

order = 2

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

trian = get_triangulation(model)
degree = order
quad = CellQuadrature(trian,degree)

strian = SkeletonTriangulation(model)
sdegree = order
squad = CellQuadrature(strian,sdegree)

V = TestFESpace(model=model,order=order,reffe=:Lagrangian,conformity=:H1,valuetype=Float64)
Q = TestFESpace(model=model,order=order-1,reffe=:Lagrangian,conformity=:L2,valuetype=Float64)

U = TrialFESpace(V)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

function a(y,x)
  v,q = y
  u,p = x
  v*u + v*p - q*p
end

function l(y)
  v,q = y
  v*4 + q
end

function a_Γ(y,x)
  v,q = y
  u,p = x
  jump(v)*mean(u) + jump(∇(q))*jump(∇(p)) - mean(v)*mean(p)
end

t_Ω = AffineFETerm(a,l,trian,quad)
t_Γ = LinearFETerm(a_Γ,strian,squad)

op = AffineFEOperator(Y,X,t_Ω,t_Γ)

op = AffineFEOperator(SparseMatrixCSR,Y,X,t_Ω,t_Γ)

op = FEOperator(Y,X,t_Ω,t_Γ)
xh = zero(X)
b = residual(op,xh)
A = jacobian(op,xh)
test_fe_operator(op,get_free_values(xh),b)

op = FEOperator(SparseMatrixCSR,Y,X,t_Ω,t_Γ)

end # module
