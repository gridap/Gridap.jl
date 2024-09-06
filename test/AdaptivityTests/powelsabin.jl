using Test
using Gridap
using Gridap.Adaptivity, Gridap.Geometry, Gridap.ReferenceFEs, Gridap.FESpaces, Gridap.Arrays
using Gridap.ReferenceFEs

using FillArrays

function l2_norm(uh,dΩ)
  sqrt(sum(∫(uh⋅uh)dΩ))
end

function l2_error(uh,u_exact,dΩ)
  eh = uh - u_exact
  return sqrt(sum(∫(eh⊙eh)dΩ))
end

Dc = 2
order = Dc

u_sol(x) = (Dc == 2) ? VectorValue(x[1],-x[2]) : VectorValue(x[1],-x[2],0.0)

domain = (Dc == 2) ? (0,1,0,1) : (0,1,0,1,0,1)
nc = (Dc == 2) ? (2,2) : (1,1,1)
model = simplexify(CartesianDiscreteModel(domain,nc))

poly  = (Dc == 2) ? TRI : TET
rrule = Adaptivity.PowellSabinRefinementRule(poly)

reffe_ref = LagrangianRefFE(Float64,poly,order)

subreffes = Fill(LagrangianRefFE(VectorValue{Dc,Float64},poly,order),Adaptivity.num_subcells(rrule))
reffe = Adaptivity.MacroReferenceFE(rrule,subreffes_u)

V = FESpace(model,reffe)

qdegree = 2*order
quad  = Quadrature(poly,Adaptivity.CompositeQuadrature(),rrule,qdegree)
Ω = Triangulation(model)
dΩ = Measure(Ω,quad)

u1(x) = VectorValue(x[1]^2 + 5*x[2]^2, x[1]^2 - 3*x[2]^2)
uh = interpolate(u1, V)
l2_error(uh,u1,dΩ)
l2_error(∇(uh),∇(u1),dΩ)

a(u,v) = ∫(u⋅v)dΩ
l(v) = ∫(u1⋅v)dΩ
op = AffineFEOperator(a,l,V,V)
uh = solve(op)
l2_error(uh,u1,dΩ)
l2_error(∇(uh),∇(u1),dΩ)
