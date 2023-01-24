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

#parent = simplexify(CartesianDiscreteModel(domain,Tuple(fill(4,D))))
#model = refine(parent;cells_to_refine=[5,12])

parent = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,Tuple(fill(2,D))))
model = refine(parent;cells_to_refine=[1,3])

trian = Triangulation(model)
ctrian = Triangulation(parent)

qorder = order*2+1
dΩ_c = Measure(ctrian,qorder)
dΩ_f = Measure(trian,qorder)

glue = get_adaptivity_glue(model)
rrules = Adaptivity.get_old_cell_refinement_rules(glue)

cell_reffe_lag = lazy_map(rr->ReferenceFE(get_polytope(rr),rr,lagrangian,Float64,order),rrules)
map(test_reference_fe,cell_reffe_lag.values)
cell_reffe_ned = lazy_map(rr->ReferenceFE(get_polytope(rr),rr,nedelec,Float64,order),rrules)
map(test_reference_fe,cell_reffe_ned.values)

reffe  = ReferenceFE(lagrangian,Float64,order)
V_c    = TestFESpace(parent,rrules,reffe;conformity=:H1,dirichlet_tags="boundary")
U_c    = TrialFESpace(V_c,sol)

V_f    = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags="boundary")
U_f    = TrialFESpace(V_f,sol)

V_c2    = TestFESpace(parent,reffe;conformity=:H1,dirichlet_tags="boundary")
U_c2    = TrialFESpace(V_c2,sol)

test_fe_space(U_c)
test_fe_space(U_f)
test_fe_space(U_c2)

# FineToCoarse interpolation, efficient due to FineToCoarseDofBasis.
u_c  = interpolate_everywhere(sol,U_c2)
u_f = interpolate_everywhere(sol,U_f)

u_f2 = Gridap.CellData.change_domain(u_c,trian,ReferenceDomain())
u_c2 = Gridap.CellData.change_domain(u_f,ctrian,ReferenceDomain())

u_fc = interpolate(u_f,U_c2)
u_fc2 = interpolate_everywhere(u_f,U_c2)

eh = u_c - u_f
e = sum(∫(eh⋅eh)*dΩ_c)

eh2 = u_c - u_fc 
e2 = sum(∫(eh2⋅eh2)*dΩ_c)

af(u,v) = ∫(v⋅u)*dΩ_f
lf(v)   = ∫(v⋅u_c)*dΩ_f
opf     = AffineFEOperator(af,lf,U_f,V_f)
u_f_pr  = solve(opf)

eh3 = u_f_pr - u_c
e3 = sum(∫(eh3⋅eh3)*dΩ_c)

sum(∫(u_f_pr)*dΩ_f)
sum(∫(u_f)*dΩ_f)
sum(∫(u_fc)*dΩ_c)


Z = collect(0.0:0.05:1.0)
for x in Z
  for y in Z
    X = VectorValue(x,y)
    print(X)
    Y = u_f(X)
    println(" -> ",Y, " - " , Y ≈ sol(X))
  end
end

M = model.model
G = M.grid

cell_coords = map(ids->G.node_coordinates[ids],G.cell_node_ids)

cv_ref = FESpaces._cell_vals(U_c,sol)
cv = FESpaces._cell_vals(U_c,u_f)

I = findall(map((a,b)->!(a≈b),cv_ref,cv))
cv[I]

i = first(I)
rr = rrules[2]

end