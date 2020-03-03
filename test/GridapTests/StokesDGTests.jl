module StokesDGTests

using Test
using Gridap
import Gridap: ∇
import LinearAlgebra: tr

const T = VectorValue{2,Float64}

u(x) = VectorValue(x[1]*x[1], x[2])
∇u(x) = TensorValue(2*x[1],0.0,0.0,1.0)
Δu(x) = VectorValue(2.0,0.0)

p(x) = x[1] - x[2]
∇p(x) = VectorValue(1.0,-1.0)

f(x) = - Δu(x) + ∇p(x)
g(x) = tr(∇u(x))

∇(::typeof(u)) = ∇u
∇(::typeof(p)) = ∇p

L = 1.0
domain = (0.0, L, 0.0, L)
ncellx = 6
partition = (ncellx,ncellx)
model = CartesianDiscreteModel(domain,partition)
model = simplexify(model)

order = 2
const h = L / ncellx
const γ = order*(order+1)
const γ0 = 1.0/10.0

V = TestFESpace(
 model=model,
 order=order,
 reffe=:Lagrangian,
 valuetype=VectorValue{2,Float64},
 conformity=:L2)

Q = TestFESpace(
 model=model,
 order=order,
 reffe=:Lagrangian,
 valuetype=Float64,
 conformity=:L2,
 constraint=:zeromean)

U = TrialFESpace(V)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

trian = get_triangulation(model)
degree = 2*order
quad = CellQuadrature(trian,degree)

btrian = BoundaryTriangulation(model)
bdegree = 2*order
bquad = CellQuadrature(btrian,bdegree)
const nb = get_normal_vector(btrian)

strian = SkeletonTriangulation(model)
sdegree = 2*order
squad = CellQuadrature(strian,sdegree)
const ns = get_normal_vector(strian)

function A_Ω(x,y)
  u, p = x
  v, q = y
  inner(∇(v), ∇(u)) - ∇(q)*u + v*∇(p)
end

function B_Ω(y)
  v, q = y
  v*f + q*g
end

function A_∂Ω(x,y)
  u, p = x
  v, q = y
  (γ/h)*v*u - v*(nb*∇(u)) - (nb*∇(v))*u + 2*(q*nb)*u
end

function B_∂Ω(y)
  v, q = y
  (γ/h)*v*u - (nb*∇(v))*u + (q*nb)*u
end

function A_Γ(x,y)
  u, p = x
  v, q = y
  (γ/h)*inner( jump(outer(v,ns)), jump(outer(u,ns))) -
    inner( jump(outer(v,ns)), mean(∇(u)) ) -
    inner( mean(∇(v)), jump(outer(u,ns)) ) +
    (γ0*h)*jump(q*ns)*jump(p*ns) +
    jump(q*ns)*mean(u) -
    mean(v)*jump(p*ns)
end

t_Ω = AffineFETerm(A_Ω,B_Ω,trian,quad)
t_∂Ω = AffineFETerm(A_∂Ω,B_∂Ω,btrian,bquad)
t_Γ = LinearFETerm(A_Γ,strian,squad)

op = AffineFEOperator(X,Y,t_Ω,t_∂Ω,t_Γ)

uh, ph = solve(op)

eu = u - uh
ep = p - ph

l2(v) = v*v
h1(v) = v*v + inner(∇(v),∇(v))

eu_l2 = sqrt(sum(integrate(l2(eu),trian,quad)))
eu_h1 = sqrt(sum(integrate(h1(eu),trian,quad)))
ep_l2 = sqrt(sum(integrate(l2(ep),trian,quad)))

tol = 1.0e-9
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

end # module
