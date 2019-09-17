module DGMultiFEOperatorsTests

# Testing Stokes problem with Dirichlet BCs for the velocity
# and pressure prescribed in a single point

using Test
using Gridap
import Gridap: ∇
import LinearAlgebra: tr

const T = VectorValue{2,Float64}

#u(x) = x[1]*x[2]*(x[1]-1)*(x[2]-1)*VectorValue(1.0,1.0)
u(x) = VectorValue(x[1]*x[1], x[2])


#function ∇u(x)
#  a = x[2]*(x[2]-1)*(2*x[1]-1)
#  b = x[1]*(x[1]-1)*(2*x[2]-1)
#  TensorValue(a,b,a,b)
#end
∇u(x) = TensorValue(2*x[1],0.0,0.0,1.0)

∇(::typeof(u)) = ∇u

#function Δu(x)
#  a = 2*x[2]*(x[2]-1) + (2*x[1]-1)*(2*x[2]-1)
#  b = 2*x[1]*(x[1]-1) + (2*x[2]-1)*(2*x[1]-1)
#  VectorValue(a,b)
#end

Δu(x) = VectorValue(2.0,0.0)

p(x) = x[1] - x[2]

∇p(x) = VectorValue(1.0,-1.0)

f(x) = - Δu(x)# + ∇p(x)

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

## Construct the FEspace 2
#fespace2 = DLagrangianFESpace(Float64,model,order)
#fixeddofs = [1,]
#fespace2 = ConstrainedFESpace(fespace2,fixeddofs)

# Define test and trial
Vh = TestFESpace(fespace1)
#Qh = TestFESpace(fespace2)
#Yh = [Vh, Qh]

Uh = TrialFESpace(fespace1)
#Ph = TrialFESpace(fespace2)
#Xh = [Uh, Ph]

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
  u = x
  v = y
  inner(∇(v), ∇(u))# - inner(∇(q), u) + inner(v, ∇(p))
end

function B_Ω(y)
  v = y
  inner(v,f)# + inner(q, g)
end

function A_∂Ω(y,x)
  u = x
  v = y
  (γ/h) * inner(v,u) - inner(v, ∇(u)*nb ) - inner(∇(v)*nb, u) #+ inner(q*nb, u) - inner(v, p*nb)
end

function B_∂Ω(y)
  v = y
  (γ/h) * inner(v,u) - inner(∇(v)*nb, ud) #+ inner(q*nb, u) - inner(v, p*nb)
end

function A_Γ(y,x)
  u = x
  v = y
  (γ/h) * inner( jump(outer(v,ns)), jump(outer(u,ns))) - inner( jump(outer(v,ns)), mean(∇(u)) ) - inner( mean(∇(v)), jump(outer(u,ns)) ) #+ (γ0*h) * inner( jump(q*ns), jump(p*ns)  ) + inner( jump(q*ns), mean(u) ) - inner( mean(v), jump(p*ns) )
end

t_Ω = AffineFETerm(A_Ω,B_Ω,trian,quad)

t_∂Ω = AffineFETerm(A_∂Ω,B_∂Ω,btrian,bquad)

t_Γ = LinearFETerm(A_Γ,strian,squad)

# Define the FEOperator
op = LinearFEOperator(Vh,Uh,t_Ω,t_∂Ω,t_Γ)

# Solve!
uh = solve(op)

uh = interpolate(Uh,u)
uh_Γ = restrict(uh,strian)
writevtk(strian,"strian",cellfields=["jump_uh" => jump(uh_Γ),"mean_uh"=>mean(uh_Γ)])

q = coordinates(squad)

ns_q = evaluate(ns,q)

phi = CellGeomap(strian)

x = evaluate(phi,q)

writevtk(x,"x",pointdata=["ns" => ns_q])


#uh = xh[1]
#ph = xh[2]

## Correct the pressure
#A = sum(integrate(p-ph,trian,quad))
#V = sum(integrate((x)->1.0,trian,quad))
#ph = ph + A/V

# Define exact solution and error
eu = u - uh

#ep = p - ph

writevtk(trian,"trian",cellfields=["uh"=>uh,"eu"=>eu])

# Define norms to measure the error
l2(u) = inner(u,u)
h1(u) = inner(∇(u),∇(u)) + l2(u)

# Compute errors
eul2 = sqrt(sum( integrate(l2(eu),trian,quad) ))
euh1 = sqrt(sum( integrate(h1(eu),trian,quad) ))

#epl2 = sqrt(sum( integrate(l2(ep),trian,quad) ))

@test eul2 < 1.e-8
@test euh1 < 1.e-8

#@test epl2 < 1.e-8


end # module
