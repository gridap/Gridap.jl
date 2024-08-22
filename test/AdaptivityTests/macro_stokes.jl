
using Gridap
using Gridap.Adaptivity, Gridap.Geometry, Gridap.ReferenceFEs, Gridap.FESpaces, Gridap.Arrays
using Gridap.ReferenceFEs

using FillArrays

# Barycentric working for D = 2,3
# PowellSabin working for D = 2 

# Parameters
Dc = 2
reftype = 2

u_sol(x) = (Dc == 2) ? VectorValue(x[1],-x[2]) : VectorValue(x[1],-x[2],0.0)
p_sol(x) = (x[1] - 1.0/2.0)

domain = (Dc == 2) ? (0,1,0,1) : (0,1,0,1,0,1)
nc = (Dc == 2) ? (4,4) : (2,2,2)
model = simplexify(CartesianDiscreteModel(domain,nc))

min_order = (reftype == 1) ? Dc : Dc-1
order = max(2,min_order)
poly  = (Dc == 2) ? TRI : TET
rrule = (reftype == 1) ? Adaptivity.BarycentricRefinementRule(poly) : Adaptivity.PowellSabinRefinementRule(poly)

reffes = Fill(LagrangianRefFE(VectorValue{Dc,Float64},poly,order),Adaptivity.num_subcells(rrule))
reffe_u = Adaptivity.MacroReferenceFE(rrule,reffes)
reffe_p = LagrangianRefFE(Float64,poly,order-1)

qdegree = 2*order
quad  = Quadrature(poly,Adaptivity.CompositeQuadrature(),rrule,qdegree)

V = FESpace(model,reffe_u,dirichlet_tags=["boundary"])
Q = FESpace(model,reffe_p,conformity=:L2,constraint=:zeromean)
U = TrialFESpace(V,u_sol)

Ω = Triangulation(model)
dΩ = Measure(Ω,quad)
@assert abs(sum(∫(p_sol)dΩ)) < 1.e-15
@assert abs(sum(∫(divergence(u_sol))dΩ)) < 1.e-15

X = MultiFieldFESpace([U,Q])
Y = MultiFieldFESpace([V,Q])

f(x) = -Δ(u_sol)(x) + ∇(p_sol)(x)
lap(u,v) = ∫(∇(u)⊙∇(v))dΩ
a((u,p),(v,q)) = lap(u,v) + ∫(divergence(u)*q - divergence(v)*p)dΩ
l((v,q)) = ∫(f⋅v)dΩ

op = AffineFEOperator(a,l,X,Y)
xh = solve(op)
uh, ph = xh
sum(∫(uh⋅uh)dΩ)
sum(∫(ph)dΩ)

eh_u = uh - u_sol
eh_p = ph - p_sol
sum(∫(eh_u⋅eh_u)dΩ)
sum(∫(eh_p*eh_p)dΩ)
