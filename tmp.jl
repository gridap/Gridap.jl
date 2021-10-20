#module tmp

using Gridap

model = CartesianDiscreteModel((0,1,0,1),(2,2))
Ω = Triangulation(model)
dΩ = Measure(Ω,2)

reffe = ReferenceFE(lagrangian,Float64,1)
V = TestFESpace(model,reffe)
U = TrialFESpace(V)

a(u,v) = ∫(∇(v)⋅∇(u))dΩ
l(v) = 0
# op = AffineFEOperator(a,l,U,V)

dv = get_fe_basis(V)
du = get_trial_fe_basis(U)

vh = FEFunction(V,rand(num_free_dofs(V)))
∇vh_q = lazy_map(Broadcasting(⋅),dv,dv)

# contribution = a(du,dv)
print_op_tree(∇vh_q,showid=true)
