module PatchAssemblersTests

using Test
using Gridap

using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs

function generate_model()
  model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(4,4)))
  labels = get_face_labeling(model)
  Geometry.add_tag_from_tags!(labels,"all",["interior","boundary"])
  merge!(
    labels,
    Geometry.face_labeling_from_cell_tags(get_grid_topology(model),[1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3],["A","B","C"])
  )
  Geometry.add_tag_from_tags!(labels,"D",["B","C"])
  return model
end

function test_patch_assembly(model,domain_tags,integration_tags)
  labels = get_face_labeling(model)

  D = num_cell_dims(model)
  Ω = Triangulation(ReferenceFE{D},model;tags=domain_tags)
  Γ = Triangulation(ReferenceFE{D-1},model;tags=domain_tags)
  
  ptopo = Geometry.PatchTopology(model;tags=domain_tags)
  n_patches = Geometry.num_patches(ptopo)
  @test n_patches == sum(Geometry.get_face_mask(labels,domain_tags,D))
  
  Ωp = Geometry.PatchTriangulation(model,ptopo;tags=integration_tags)
  Γp = Geometry.PatchBoundaryTriangulation(model,ptopo;tags=integration_tags)
  @test Geometry.num_cells(Ωp) == sum(Geometry.get_face_mask(labels,integration_tags,D))
  @test length(unique(get_glue(Γp,Val(1)).tface_to_mface)) == sum(Geometry.get_face_mask(labels,integration_tags,D-1))
  
  order = 1
  VΩ = FESpaces.PolytopalFESpace(Ω,Float64,order)
  VΓ = FESpaces.PolytopalFESpace(Γ,Float64,order)
  X = MultiFieldFESpace([VΩ,VΓ])
  
  qdegree = 2*(order+1)
  dΩp = Measure(Ωp,qdegree)
  dΓp = Measure(Γp,qdegree)
  
  f(x) = 1.0
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v,Ω,dΩ) = ∫(u⋅Π(v,Ω))dΩ
  laplacian(u,v,dΩ) = ∫(∇(u)⋅∇(v))dΩ
  
  passem1 = FESpaces.PatchAssembler(ptopo,VΩ,VΩ)
  a1(uΩ,vΩ) = laplacian(uΩ,vΩ,dΩp)
  l1(vΩ) = ∫(f⋅vΩ)dΩp
  matvecs1 = assemble_matrix_and_vector(a1,l1,passem1,VΩ,VΩ)
  @test length(collect(matvecs1)) == n_patches
  
  passem2 = FESpaces.PatchAssembler(ptopo,VΓ,VΓ)
  a2(uΓ,vΓ) = mass(uΓ,vΓ,Γp,dΓp)
  l2(vΓ) = ∫(f⋅vΓ)dΓp
  matvecs2 = assemble_matrix_and_vector(a2,l2,passem2,VΓ,VΓ)
  @test length(collect(matvecs2)) == n_patches
  
  passem3 = FESpaces.PatchAssembler(ptopo,X,X)
  a3((uΩ,uΓ),(vΩ,vΓ)) = laplacian(uΩ,vΩ,dΩp) + mass(uΓ,vΩ,Γp,dΓp) + mass(uΩ,vΓ,Γp,dΓp) + mass(uΓ,vΓ,Γp,dΓp)
  l3((vΩ,vΓ)) = ∫(f⋅vΩ)dΩp
  matvecs3 = assemble_matrix_and_vector(a3,l3,passem3,X,X)
  @test length(collect(matvecs3)) == n_patches
end

model = generate_model()

# Case 0: 
#  - PatchTopology spans all cells
#  - Integration on all patches
test_patch_assembly(model,"all","all")

# Case 1:
#  - PatchTopology spans all cells
#  - Integration a subset of patches
test_patch_assembly(model,"all","D")

# Case 2:
#  - PatchTopology spans a subset of cells
#  - Integration on all patches
test_patch_assembly(model,"D","D")

# Case 3:
#  - PatchTopology spans a subset of cells
#  - Integration a subset of patches
test_patch_assembly(model,"D","B")

end