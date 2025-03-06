
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.Arrays

using DataStructures

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
x1 = lazy_map(FESpaces.LocalSolveMap(),A,b)
x2 = lazy_map(FESpaces.LocalSolveMap(),Ab)

# MultiField + Static condensation

using Gridap.MultiField
using Gridap.CellData

Π(u) = change_domain(u,Γp,DomainStyle(u))
function a_mf((u1,_u2),(v1,_v2))
  u2, v2 = Π(_u2), Π(_v2)
  aΩp(u1,v1) + aΓp(u2,v2) + aΓp(u1,v2) + aΓp(u2,v1)
end
function l_mf((v1,_v2))
  v2 = Π(_v2)
  lΩp(v1) + 2*lΓp(v2) + 2*lΓp(v1)
end

mfs = MultiField.ConsecutiveMultiFieldStyle()
X = MultiFieldFESpace([VΩ,VΓ];style=mfs)
assem = FESpaces.PatchAssembler(ptopo,X,X)
Ab = assemble_matrix_and_vector(a_mf,l_mf,assem,X,X)
x = lazy_map(FESpaces.LocalSolveMap(),Ab)

mfs = MultiField.BlockMultiFieldStyle()
X = MultiFieldFESpace([VΩ,VΓ];style=mfs)
assem = FESpaces.PatchAssembler(ptopo,X,X)
Ab = assemble_matrix_and_vector(a_mf,l_mf,assem,X,X)
x = lazy_map(FESpaces.StaticCondensationMap(),Ab)

########################################################################

model = CartesianDiscreteModel((0,1,0,1),(4,4))

ptopo = Geometry.PatchTopology(model)

patch_cells = Geometry.get_patch_cells(ptopo)
patch_facets = Geometry.get_patch_facets(ptopo)
patch_nodes = Geometry.get_patch_faces(ptopo,0)

Ωp = Geometry.PatchTriangulation(model,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)
Λp = Geometry.PatchTriangulation(ReferenceFE{1},model,ptopo)

Ω = Triangulation(model)
Γ = Triangulation(ReferenceFE{1},model)

order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
VΩ = FESpace(Ω,reffe,conformity=:L2)
VΓ = FESpace(Γ,reffe,conformity=:L2)

mface_to_dofs = get_cell_dof_ids(VΓ,Triangulation(ReferenceFE{1},model))
pface_to_dofs = get_cell_dof_ids(VΓ,Λp)
patch_to_dofs = Arrays.merge_entries(mface_to_dofs,patch_to_lpface_to_mface;acc=OrderedSet{Int32}())
pface_to_ldofs = find_local_index(pface_to_dofs,pface_to_patch,patch_to_dofs)

Df = 1
patch_to_lpface_to_mface = Geometry.get_patch_faces(ptopo,Df)
pface_to_patch = Geometry.get_pface_to_patch(ptopo,Df)
pface_to_lpface = Geometry.get_pface_to_lpface(ptopo,Df)

strian = Γp # get_triangulation(VΓ)
sglue = get_glue(strian,Val(Df))
sface_to_mface = sglue.tface_to_mface
mface_to_sface = sglue.mface_to_tface

sface_to_pface = strian.tface_to_pface
sface_to_patch = lazy_map(Reindex(pface_to_patch),sface_to_pface)
sface_to_lpface = lazy_map(Reindex(pface_to_lpface),sface_to_pface)


V = FESpace(Γ,reffe,conformity=:L2,dirichlet_tags="boundary")

nfree = num_free_dofs(V)
ndir = num_dirichlet_dofs(V)
free_ids = ndir+1:nfree+ndir
dir_ids = ndir:-1:1

V0 = FESpaces.renumber_free_and_dirichlet_dof_ids(V,free_ids,dir_ids)

v_ids = get_cell_dof_ids(V)
v0_ids = get_cell_dof_ids(V0)

p_v_ids = FESpaces.generate_patch_dof_ids(V,ptopo)
p_v0_ids = FESpaces.generate_patch_dof_ids(V0,ptopo)

pface_to_pdof_v = FESpaces.generate_pface_to_pdofs(V,ptopo,Λp,p_v_ids)
pface_to_pdof_v0 = FESpaces.generate_pface_to_pdofs(V0,ptopo,Λp,p_v0_ids)
pface_to_pdof_v == pface_to_pdof_v0


PV = FESpaces.PatchFESpace(V,ptopo)


