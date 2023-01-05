module FineToCoarseFieldsTests

using Test
using Gridap
using Gridap.Algebra
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using Gridap.FESpaces
using FillArrays

D = 2
order  = 1
qorder = order*2+1
domain = Tuple(repeat([0,1],D))

parent = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,Tuple(fill(4,D))))
model = refine(parent;cells_to_refine=[1,2])

trian = Triangulation(model)
ctrian = Triangulation(parent)

glue = get_adaptivity_glue(model)
rrules = Adaptivity.get_old_cell_refinement_rules(glue)

cell_quad = lazy_map(rr -> Quadrature(rr,qorder;bundle_points=true),rrules)

dΩ_comp = Measure(ctrian,rrules,qorder;bundle_points=true)

reffe = ReferenceFE(lagrangian,Float64,order)
V_f = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags="boundary")
U_f = TrialFESpace(V_f,sol)
V_c = TestFESpace(parent,reffe;conformity=:H1,dirichlet_tags="boundary")
U_c = TrialFESpace(V_c,sol)

# Fine FEFunction -> Coarse FEFunction, by projection
ac(u,v) = ∫(v⋅u)*dΩ_c
lc(v)   = ∫(v⋅uh_f_inter)*dΩ_comp
opc     = AffineFEOperator(ac,lc,U_c,V_c)
uh_c_pr = solve(opc)

v_c_pr = map(p -> uh_c_pr(p), pts)
@test v_c_pr ≈ v_r

eh = sum(∫(uh_f_inter-uh_c_pr)*dΩ_c)
@test eh < 1.e8


end