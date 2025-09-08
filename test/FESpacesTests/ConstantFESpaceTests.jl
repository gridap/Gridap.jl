module ConstantFESpacesTests

using Gridap
using Gridap.FESpaces
using Test

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)
Λ = ConstantFESpace(model)
Gridap.FESpaces.test_fe_space(Λ)
M = TrialFESpace(Λ)

order = 2
u((x,y)) = (x+y)^order
f(x) = -Δ(u,x)
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe,dirichlet_tags="boundary")
U = TrialFESpace(u,V)
Y = MultiFieldFESpace([V,Λ])
X = MultiFieldFESpace([U,M])

Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)
a((u,μ),(v,λ)) = ∫(∇(v)⋅∇(u))dΩ + ∫(v*μ)dΩ + ∫(λ*u)dΩ
l((v,λ)) = ∫(λ*u)dΩ + ∫(v*f)dΩ
op = AffineFEOperator(a,l,X,Y)
uh = solve(op)

@assert sum(∫((uh[1]-u)*(uh[1]-u))dΩ) < 1.0e-14
abs(sum(∫(uh[2])dΩ)) < 1.0e-12

Λ2 = ConstantFESpace(model,field_type=VectorValue{2,Float64})
Gridap.FESpaces.test_fe_space(Λ2)
M2 = TrialFESpace(Λ2)
a2(μ,λ) = ∫(λ⋅μ)dΩ
l2(λ) = ∫(VectorValue(0.0,0.0)⋅λ)dΩ
op2 = AffineFEOperator(a2,l2,M2,Λ2)
μ2h = solve(op2)
@assert sum(∫(μ2h⋅μ2h)dΩ) < 1.0e-12

trian = Triangulation(model,[1,2,3,4])
Λ3 = ConstantFESpace(trian,field_type=VectorValue{2,Float64})
Gridap.FESpaces.test_fe_space(Λ3)

# MultiConstantFESpace

tags = ["tag_5","tag_6","tag_7","tag_8"]
Λ4 = FESpaces.MultiConstantFESpace(model,tags,1)
FESpaces.test_fe_space(Λ4)

btrians = map(tag -> Boundary(model,tags=tag), tags)
Λ5 = FESpaces.MultiConstantFESpace(btrians)
FESpaces.test_fe_space(Λ5)

t4 = get_triangulation(Λ4)
t5 = get_triangulation(Λ5)
@test num_cells(t4) == num_cells(t5)

function sorted_dof_ids(space)
  trian = get_triangulation(space)
  tface_to_mface = Gridap.Geometry.get_glue(trian,Val(num_cell_dims(trian))).tface_to_mface
  cell_dof_ids = get_cell_dof_ids(space)
  return cell_dof_ids[sortperm(tface_to_mface)]
end

@test sorted_dof_ids(Λ4) == sorted_dof_ids(Λ5)

end # module
