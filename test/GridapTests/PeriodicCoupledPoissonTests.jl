module PeriodicCoupledPoissonTests

using Gridap
using Test
using LinearAlgebra

u(x) =  x[1]^2 + 2*x[2]^2
v(x) = -x[2]^2
f(x) = -Δ(u)(x)
g(x) = -Δ(v)(x) - u(x)

model = CartesianDiscreteModel((0,1,0,1,0,0.01),(4,4,3);isperiodic=(false,false,true))
order = 2
labels = get_face_labeling(model)

add_tag_from_tags!(labels,"dirichlet",append!(collect(1:20),[23,24,25,26]))

trian = get_triangulation(model)
degree = 2*order
quad = CellQuadrature(trian,degree)

Vu = FESpace(
       reffe=:Lagrangian, order=order, valuetype=Float64,
       conformity=:H1, model=model, dirichlet_tags="dirichlet")
Vv = FESpace(
      reffe=:Lagrangian, order=order, valuetype=Float64,
      conformity=:H1, model=model, dirichlet_tags="dirichlet")

U = TrialFESpace(Vu,u)
V = TrialFESpace(Vv,v)

X = MultiFieldFESpace([U, V])
Y = MultiFieldFESpace([Vu, Vv])

function aa(X,Y)
  u, v = X
  v_u, v_v = Y
  inner(∇(v_u),∇(u)) + inner(∇(v_v),∇(v)) - v_v*u
end

function l(Y)
  v_u, v_v = Y
  v_u*f + v_v*g
end

t_Ω = AffineFETerm(aa,l,trian,quad)
op = AffineFEOperator(X,Y,t_Ω)
uh, vh = solve(op)


eu = u - uh
ev = v - vh

l2(u) = u*u
h1(u) = ∇(u)⋅∇(u) + l2(u)

eul2 = sqrt(sum( integrate(l2(eu),trian,quad) ))
euh1 = sqrt(sum( integrate(h1(eu),trian,quad) ))
ul2 = sqrt(sum( integrate(l2(uh),trian,quad) ))
uh1 = sqrt(sum( integrate(h1(uh),trian,quad) ))

evl2 = sqrt(sum( integrate(l2(ev),trian,quad) ))
evh1 = sqrt(sum( integrate(h1(ev),trian,quad) ))
vl2 = sqrt(sum( integrate(l2(vh),trian,quad) ))
vh1 = sqrt(sum( integrate(h1(vh),trian,quad) ))

@test eul2/ul2 < 1.e-8
@test euh1/uh1 < 1.e-7
@test evl2/vl2 < 1.e-8
@test evh1/vh1 < 1.e-7


end #module
