module FineToCoarseFieldsTests

using Test
using Gridap
using Gridap.Arrays
using Gridap.Algebra
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using Gridap.FESpaces
using FillArrays

sol(x) = x[1] + x[2]

D = 2
order  = 1
qorder = order*2+1
domain = Tuple(repeat([0,1],D))

parent = CartesianDiscreteModel(domain,Tuple(fill(2,D)))
model  = refine(parent)

trian  = Triangulation(model)
ctrian = Triangulation(parent)

qorder = order*2+1
dΩ_c = Measure(ctrian,qorder)
dΩ_f = Measure(trian,qorder)
dΩ_cf = Measure(ctrian,trian,qorder)

glue = get_adaptivity_glue(model)
rrules = Adaptivity.get_old_cell_refinement_rules(glue)

cell_reffe_lag = lazy_map(rr->ReferenceFE(get_polytope(rr),rr,lagrangian,Float64,order),rrules)
lazy_map(test_reference_fe,cell_reffe_lag)
cell_reffe_ned = lazy_map(rr->ReferenceFE(get_polytope(rr),rr,nedelec,Float64,order),rrules)
lazy_map(test_reference_fe,cell_reffe_ned)

# Lagrangian tests
reffe  = ReferenceFE(lagrangian,Float64,order)
V_c    = TestFESpace(parent,rrules,reffe;conformity=:H1,dirichlet_tags="boundary")
U_c    = TrialFESpace(V_c,sol)

V_f    = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags="boundary")
U_f    = TrialFESpace(V_f,sol)

test_fe_space(U_c)
test_fe_space(U_f)

u_c = interpolate_everywhere(sol,U_c)
u_f = interpolate_everywhere(sol,U_f)

u_fc = interpolate(u_f,U_c)
u_fc2 = interpolate_everywhere(u_f,U_c)

eh = u_c - u_f
@test sum(∫(eh⋅eh)*dΩ_f) < 1.e-12

eh2 = u_c - u_fc
@test sum(∫(eh2⋅eh2)*dΩ_c) < 1.e-12

eh3 = u_c - u_fc2 
@test sum(∫(eh3⋅eh3)*dΩ_c) < 1.e-12

# Moment-based (Nedelec) tests
sol((x,y)) = 2*VectorValue(-y,x)

reffe  = ReferenceFE(nedelec,Float64,order)
V_c    = TestFESpace(parent,rrules,reffe;dirichlet_tags="boundary")
U_c    = TrialFESpace(V_c,sol)

V_f    = TestFESpace(model,reffe;dirichlet_tags="boundary")
U_f    = TrialFESpace(V_f,sol)

test_fe_space(U_c)
test_fe_space(U_f)

u_c = interpolate_everywhere(sol,U_c)
u_f = interpolate_everywhere(sol,U_f)

u_fc = interpolate(u_f,U_c)
u_fc2 = interpolate_everywhere(u_f,U_c)

eh = u_c - u_f
@test sum(∫(eh⋅eh)*dΩ_f) < 1.e-12

eh2 = u_c - u_fc
@test sum(∫(eh2⋅eh2)*dΩ_c) < 1.e-12

eh3 = u_c - u_fc2 
@test sum(∫(eh3⋅eh3)*dΩ_c) < 1.e-12

modelH=CartesianDiscreteModel((0,1,0,1),(1,1))
modelh=refine(modelH,2)
reffe=LagrangianRefFE(Float64,QUAD,1)
XH  = TestFESpace(modelH,reffe)
xH  = get_fe_basis(XH)
xHh = change_domain(xH,get_triangulation(modelh),ReferenceDomain())
evaluate(Gridap.CellData.get_data(xHh)[1],
         [Point(0.0,0.0),Point(0.5,0.5)])

end