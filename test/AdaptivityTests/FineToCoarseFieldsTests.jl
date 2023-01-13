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

parent = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,Tuple(fill(4,D))))
model = refine(parent;cells_to_refine=[1,2])

trian = Triangulation(model)
ctrian = Triangulation(parent)

glue = get_adaptivity_glue(model)
rrules = Adaptivity.get_old_cell_refinement_rules(glue)

reffe  = ReferenceFE(lagrangian,Float64,order)
V_c    = TestFESpace(parent,rrules,reffe;conformity=:H1)
U_c    = TrialFESpace(V_c,sol)

V_f    = TestFESpace(model,reffe;conformity=:H1)
U_f    = TrialFESpace(V_f,sol)

# FineToCoarse interpolation, efficient due to FineToCoarseDofBasis. 
u_f  = interpolate(sol,U_f)
u_fc = interpolate(u_f,U_c)



end