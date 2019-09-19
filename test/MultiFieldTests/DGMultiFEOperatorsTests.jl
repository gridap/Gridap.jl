module DGMultiFEOperatorsTests

# Testing Stokes problem with Dirichlet BCs for the velocity
# and pressure prescribed in a single point

using Test
using Gridap
import Gridap: ∇
import LinearAlgebra: tr

const T = VectorValue{2,Float64}

u(x) = VectorValue(x[1]*x[1], x[2])

∇u(x) = TensorValue(2*x[1],0.0,0.0,1.0)

∇(::typeof(u)) = ∇u

Δu(x) = VectorValue(2.0,0.0)

p(x) = x[1] - x[2]

∇p(x) = VectorValue(1.0,-1.0)

f(x) = - Δu(x) + ∇p(x)

g(x) = tr(∇u(x))

ud(x) = u(x)

# Discrete model
L = 1.0
limits = (0.0, L, 0.0, L)
ncellx = 4
model = CartesianDiscreteModel(domain=limits, partition=(ncellx,ncellx))
model = simplexify(model)

h = L / ncellx

order = 2
γ = order*(order+1)
γ0 = 1.0/10.0

# Construct the FEspace 1
fespace1 = DLagrangianFESpace(T,model,order)

# Construct the FEspace 2
fespace2 = DLagrangianFESpace(Float64,model,order)
fespace2 = ZeroMeanFESpace(fespace2,order)

# Define test and trial
Vh = TestFESpace(fespace1)
Qh = TestFESpace(fespace2)
Yh = [Vh, Qh]

Uh = TrialFESpace(fespace1)
Ph = TrialFESpace(fespace2)
Xh = [Uh, Ph]

# Define integration mesh and quadrature for volume
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2*order)

btrian = BoundaryTriangulation(model)
bquad = CellQuadrature(btrian,order=2*order)
nb = NormalVector(btrian)

strian = SkeletonTriangulation(model)
squad = CellQuadrature(strian,order=2*order)
ns = NormalVector(strian)

# Weak form

function A_Ω(y,x)
  u, p = x
  v, q = y
  inner(∇(v), ∇(u)) - inner(∇(q), u) + inner(v, ∇(p))
end

function B_Ω(y)
  v, q = y
  inner(v,f) + inner(q, g)
end

function A_∂Ω(y,x)
  u, p = x
  v, q = y
  (γ/h) * inner(v,u) - inner(v, ∇(u)*nb )  - inner(∇(v)*nb, u) + 2*inner(q*nb,u) # - inner(v, p*nb) + inner(v,p*nb)
end

function B_∂Ω(y)
  v, q = y
  (γ/h) * inner(v,ud) - inner(∇(v)*nb, ud) + inner(q*nb, ud)
end

function A_Γ(y,x)
  u, p = x
  v, q = y
  (γ/h) * inner( jump(outer(v,ns)), jump(outer(u,ns))) -
    inner( jump(outer(v,ns)), mean(∇(u)) ) -
    inner( mean(∇(v)), jump(outer(u,ns)) ) +
    (γ0*h) * inner( jump(q*ns), jump(p*ns)  ) +
    inner( jump(q*ns), mean(u) ) -
    inner( mean(v), jump(p*ns) )
end

t_Ω = AffineFETerm(A_Ω,B_Ω,trian,quad)

t_∂Ω = AffineFETerm(A_∂Ω,B_∂Ω,btrian,bquad)

t_Γ = LinearFETerm(A_Γ,strian,squad)

# Define the FEOperator
op = LinearFEOperator(Yh,Xh,t_Ω,t_∂Ω,t_Γ)

# Solve!
xh = solve(op)
uh, ph = xh

# Define exact solution and error
eu = u - uh

ep = p - ph

#writevtk(trian,"trian",cellfields=["uh"=>uh,"ph"=>ph, "eu"=>eu, "ep"=>ep])

# Define norms to measure the error
l2(u) = inner(u,u)
h1(u) = inner(∇(u),∇(u)) + l2(u)

# Compute errors
eul2 = sqrt(sum( integrate(l2(eu),trian,quad) ))
euh1 = sqrt(sum( integrate(h1(eu),trian,quad) ))

epl2 = sqrt(sum( integrate(l2(ep),trian,quad) ))

@test eul2 < 1.e-8
@test euh1 < 1.e-8

@test epl2 < 1.e-8


end # module
