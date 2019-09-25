module DarcyFEOperatorsTests

##
using Test
using Gridap
import Gridap: ∇

u(x) = VectorValue(x[1]*x[2], -0.5*x[2]^2)
∇u(x) = TensorValue(x[2],0.0,x[1],-x[2])
divu(x) = 0.0
∇(::typeof(u)) = ∇u

p(x) = x[1] + x[2]
∇p(x) = VectorValue(1.0,1.0)
∇(::typeof(p)) = ∇p

f(x) = divu(x)
g(x) = p(x)

const kinv_1 = TensorValue(1.0,0.0,0.0,1.0)
r(x) = kinv_1*u(x) + ∇p(x)

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(50,50))
order = 2
V = FESpace( reffe=:RaviartThomas, conformity=:HDiv, order=2, model=model, diritags = [5,6])
_Q = FESpace( reffe=:QLagrangian, conformity=:L2, valuetype = Float64, order = 1,
                    model = model)
# If all faces of the domain are zero values, we must fix zero mean value
# _Q = FESpace( reffe=:QLagrangian, conformity=:L2, valuetype = Float64, order = 1,
                    # model = model, constraint = :zeromean)

V_0 = TestFESpace(V)
Q = TestFESpace(_Q)
Y = [V_0, Q]

V_g = TrialFESpace(V,u)
X = [V_g, Q]

trian = Triangulation(model)
quad = CellQuadrature(trian,degree=2)

@law σ(x,u) = kinv_1*u

function a(y,x)
  v, q = y
  u, p = x
  inner(v,σ(u)) - inner(div(v),p) + inner(q,div(u))
end

function b_Ω(y)
  v, q = y
  inner(v,r) + inner(q,f)
end

t_Ω = AffineFETerm(a,b_Ω,trian,quad)

neumanntags = [7,8]
btrian = BoundaryTriangulation(model,neumanntags)
bquad = CellQuadrature(btrian,degree=order*2)
nb = NormalVector(btrian)

function b_Γ(y)
  v, q = y
  - inner(v*nb,g)
end

t_Γ = FESource(b_Γ,btrian,bquad)

op = LinearFEOperator(Y,X,t_Ω,t_Γ)

uh, ph = solve(op)

e1 = u - uh
e2 = p - ph

# Define norms to measure the error
l2(u) = inner(u,u)
hdiv(u) = inner(div(u),div(u)) + l2(u)

# Compute errors
e1l2 = sqrt(sum( integrate(l2(e1),trian,quad) ))
e1hdiv = sqrt(sum( integrate(hdiv(e1),trian,quad) ))
e2l2 = sqrt(sum( integrate(l2(e2),trian,quad) ))

# writevtk(trian,"/home/santiago/github-repos/Gridap/tmp/darcyresults",
         # cellfields=["uh"=>uh,"ph"=>ph, "euh"=> e1, "eph"=> e2])

@test e1l2 < 1.e-8
@test e1hdiv < 1.e-8
@test e2l2 < 1.e-8
##
end # module
