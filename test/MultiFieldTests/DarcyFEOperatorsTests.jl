module DarcyFEOperatorsTests

# Testing Darvy problem with Dirichlet BCs for the velocity
# and pressure prescribed in a single point
##
using Test
using Gridap
import Gridap: ∇

const T = VectorValue{2,Float64}

# Define manufactured functions
u(x) = VectorValue(x[1], -1*x[2])
∇u(x) = TensorValue(1.0,0.0,0.0,1.0)
∇(::typeof(u)) = ∇u


# @santiagobadia :  Bug, if p(x) is equal to zero on the Neumann boundary
# everything works. When it is different from zero, it does not work
p(x) = x[1] # It works
# p(x) = x[1] - 0.5 # It does not work
∇p(x) = VectorValue(1.0,0.0)
∇(::typeof(p)) = ∇p


# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(3,3))

# Construct the FEspace 1
order = 2
fespace1 = FESpace( reffe=:RaviartThomas, conformity=:HDiv, order=2, model=model, diritags = [5,6,8])
# If all faces of the domain are zero values, we must fix zero mean value
# fespace2 = FESpace( reffe=:PLagrangian, conformity=:L2, valuetype = Float64, order = 1,
                    # model = model, constraint = :zeromean)
fespace2 = FESpace( reffe=:PLagrangian, conformity=:L2, valuetype = Float64, order = 1,
                    model = model)

# Define test and trial
V1 = TestFESpace(fespace1)
V2 = TestFESpace(fespace2)
V = [V1, V2]

U1 = TrialFESpace(fespace1,u)
U2 = TrialFESpace(fespace2)
U = [U1, U2]

# Define integration mesh and quadrature for volume
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

const kinv_1 = TensorValue(1.0,0.0,0.0,1.0)
@law σ(x,u) = kinv_1*u

function a(y,x)
  v, q = y
  u, p = x
  inner(v,σ(u)) - inner(div(v),p) + inner(q,div(u))
end

function b_Ω(y)
  v, q = y
  inner(v,u) + inner(v,∇(p))
end

t_Ω = AffineFETerm(a,b_Ω,trian,quad)

neumanntags = [7]
btrian = BoundaryTriangulation(model,neumanntags)
bquad = CellQuadrature(btrian,order=order*2)
nb = NormalVector(btrian)

function b_Γ(y)
  v, q = y
  inner(v*nb,p)
end

t_Γ = FESource(b_Γ,btrian,bquad)

# op = LinearFEOperator(V,U,t_Ω)
op = LinearFEOperator(V,U,t_Ω,t_Γ)

# Solve!
uh = solve(op)

# Define exact solution and error
e1 = u - uh[1]

e2 = p - uh[2]

pi = interpolate(fespace2,p)
m(p) = inner(p,1.0)
sum(integrate(m(pi),trian,quad))

# Define norms to measure the error
l2(u) = inner(u,u)
# h1(u) = inner(∇(u),∇(u)) + l2(u)

# Compute errors
e1l2 = sqrt(sum( integrate(l2(e1),trian,quad) ))
# e1h1 = sqrt(sum( integrate(h1(e1),trian,quad) ))

e2l2 = sqrt(sum( integrate(l2(e2),trian,quad) ))

@test e1l2 < 1.e-8
# @test e1h1 < 1.e-8

@test e2l2 < 1.e-8
# writevtk(trian,"/home/santiago/github-repos/Gridap/tmp/darcyresults",
        # cellfields=["uh"=>uh[1],"ph"=>uh[2], "pe"=>pi, "eh"=> e2])
##
end # module
