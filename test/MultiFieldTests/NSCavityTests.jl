using Gridap
using Test
n = 10
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(n,n))
labels = FaceLabels(model)
add_tag_from_tags!(labels,"diri1",[6,])
add_tag_from_tags!(labels,"diri0",[1,2,3,4,5,7,8])
D = 2
order = 2
fespace1 = FESpace(
  reffe=:Lagrangian, conformity=:H1, valuetype=VectorValue{D,Float64},
  model=model, labels=labels, order=order, diritags=["diri0","diri1"])
fespace2 = FESpace(
  reffe=:PLagrangian, conformity=:L2, valuetype=Float64,
  model=model, order=order-1, constraint=:zeromean)
uD0(x) = VectorValue(0.0,0.0)
uD1(x) = VectorValue(1.0,0.0)

V = TestFESpace(fespace1)
Q = TestFESpace(fespace2)
Y = [V, Q]

U = TrialFESpace(fespace1,[uD0,uD1])
P = TrialFESpace(fespace2)
X = [U, P]
const Re = 10.0
@law conv(x,u,∇u) = Re*(∇u')*u
@law dconv(x,du,∇du,u,∇u) = conv(x,u,∇du)+conv(x,du,∇u)

function a(y,x)
  u, p = x
  v, q = y
  inner(∇(v),∇(u)) - inner(div(v),p) + inner(q,div(u))
end

c(v,u) = inner(v,conv(u,∇(u)))
dc(v,du,u) = inner(v,dconv(du,∇(du),u,∇(u)))

function res(x,y)
  u, p = x
  v, q = y
  a(y,x) + c(v,u)
end

function jac(x,y,dx)
  u, p = x
  v, q = y
  du, dp = dx
  a(y,dx)+ dc(v,du,u)
end
trian = Triangulation(model)
quad = CellQuadrature(trian,degree=(order-1)*2)
t_Ω = NonLinearFETerm(res,jac,trian,quad)
op = NonLinearFEOperator(Y,X,t_Ω)
using LineSearches: BackTracking
nls = NLSolver(
  show_trace=false, method=:newton, linesearch=BackTracking())
solver = NonLinearFESolver(nls)
uh, ph = solve(solver,op)
#writevtk(trian,"ins-results",cellfields=["uh"=>uh,"ph"=>ph])

meanval = sum( integrate(ph,trian,quad))
@test meanval + 1 ≈ 1

