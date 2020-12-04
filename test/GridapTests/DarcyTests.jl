module DarcyTests

using Test
using Gridap
import Gridap: ∇, divergence
using LinearAlgebra

u(x) = VectorValue(2*x[1],x[1]+x[2])

divergence(::typeof(u)) = (x) -> 3

p(x) = x[1]-x[2]

∇p(x) = VectorValue(1,-1)

∇(::typeof(p)) = ∇p

f(x) = u(x) + ∇p(x)

domain = (0,1,0,1)
partition = (4,4)
order = 1
model = CartesianDiscreteModel(domain,partition)

V = FESpace(
  reffe=:RaviartThomas, order=order, valuetype=VectorValue{2,Float64},
  conformity=:Hdiv, model=model, dirichlet_tags=[5,6])

Q = FESpace(
  reffe=:QLagrangian, order=order, valuetype=Float64,
  conformity=:L2, model=model)

U = TrialFESpace(V,u)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

trian = Triangulation(model)
degree = 2
quad = CellQuadrature(trian,degree)

neumanntags = [7,8]
btrian = BoundaryTriangulation(model,tags=neumanntags)
degree = 2*(order+1)
bquad = CellQuadrature(btrian,degree)
nb = get_normal_vector(btrian)

function a(x,y)
  u, p = x
  v, q = y
  u⋅v - p*(∇⋅v) + q*(∇⋅u)
end

function l(y)
  v, q = y
  v⋅f + q*(∇⋅u)
end

function l_Γ(y)
  v, q = y
  -(v⋅nb)*p
end

t_Ω = AffineFETerm(a,l,trian,quad)
t_Γ = FESource(l_Γ,btrian,bquad)
op = AffineFEOperator(X,Y,t_Ω,t_Γ)
xh = solve(op)
uh, ph = xh

eu = u - uh
ep = p - ph

l2(v) = v⋅v
h1(v) = v*v + ∇(v)⋅∇(v)

eu_l2 = sum(integrate(l2(eu),trian,quad))
ep_l2 = sum(integrate(l2(ep),trian,quad))
ep_h1 = sum(integrate(h1(ep),trian,quad))

tol = 1.0e-9
@test eu_l2 < tol
@test ep_l2 < tol
@test ep_h1 < tol

end # module
