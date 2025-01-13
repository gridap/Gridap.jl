
using Gridap
using Gridap.Geometry, Gridap.FESpaces

model = CartesianDiscreteModel((0,1,0,1),(4,4))

ptopo = Geometry.PatchTopology(model)

patch_cells = Geometry.get_patch_cells(ptopo)
patch_facets = Geometry.get_patch_facets(ptopo)
patch_nodes = Geometry.get_patch_faces(ptopo,0)

Ωp = Geometry.PatchTriangulation(model,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)

Ω = Triangulation(model)
Γ = Triangulation(ReferenceFE{1},model)

order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
VΩ = FESpace(Ω,reffe)
VΓ = FESpace(Γ,reffe)

assem_Ω = FESpaces.PatchAssembler(ptopo,VΩ,VΩ)
assem_Γ = FESpaces.PatchAssembler(ptopo,VΓ,VΓ)

f(x) = x[1] + x[2]
biform(u,v,dΩ) = ∫(u⋅v)dΩ
liform(v,dΩ) = ∫(v⋅f)dΩ
aΩ(u,v) = biform(u,v,Measure(Ω,2*order))
aΓ(u,v) = biform(u,v,Measure(Γ,2*order))

aΩp(u,v) = biform(u,v,Measure(Ωp,2*order))
aΓp(u,v) = biform(u,v,Measure(Γp,2*order))
a_mixed(u,v) = aΩp(u,v) + aΓp(u,v)

lΩp(v) = liform(v,Measure(Ωp,2*order))
lΓp(v) = liform(v,Measure(Γp,2*order))
l_mixed(v) = lΩp(v) + lΓp(v)

AΩ = assemble_matrix(aΩp,assem_Ω,VΩ,VΩ)
bΩ = assemble_vector(lΩp,assem_Ω,VΩ)
AbΩ = assemble_matrix_and_vector(aΩp,lΩp,assem_Ω,VΩ,VΩ)
AΩ_arr = collect_cell_matrix(VΩ,VΩ,aΩ(get_fe_basis(VΩ),get_trial_fe_basis(VΩ)))[1][1]
all(map(isequal,AΩ,AΩ_arr))

AΓ = assemble_matrix(aΓp,assem_Γ,VΓ,VΓ)
AΓ_arr = collect_cell_matrix(VΓ,VΓ,aΓ(get_fe_basis(VΓ),get_trial_fe_basis(VΓ)))[1][1]

A = assemble_matrix(a_mixed,assem_Ω,VΩ,VΩ)
b = assemble_vector(l_mixed,assem_Ω,VΩ)
Ab = assemble_matrix_and_vector(a_mixed,l_mixed,assem_Ω,VΩ,VΩ)
