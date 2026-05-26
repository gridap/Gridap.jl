module PatchAssemblersTests

using Test
using Gridap

using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs
using Gridap.Arrays

function projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  P = LocalOperator(
    LocalSolveMap(), V0, mass, mass; trian_out = Ω
  )
  return P
end

function reconstruction_operator(ptopo,X,Ω,Γp,dΩp,dΓp,order)
  reffe_L = ReferenceFE(lagrangian, Float64, order+1; space=:P)
  L = FESpace(Ω, reffe_L; conformity=:L2)
  reffe_Λ = ReferenceFE(lagrangian, Float64, 0; space=:P)
  Λ = FESpace(Ω, reffe_Λ; conformity=:L2)

  n = get_normal_vector(Γp)
  Πn(v) = ∇(v)⋅n
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((u,λ),(v,μ))   =  ∫( (∇(u)⋅∇(v)) + (μ*u) + (λ*v) )dΩp
  rhs((uT,uF),(v,μ)) =  ∫( (∇(uT)⋅∇(v)) + (uT*μ) )dΩp + ∫( (uF - Π(uT,Γp))*(Πn(v)) )dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  W = MultiFieldFESpace([L,Λ];style=MultiField.MultiFieldStyle(Y))
  R = LocalOperator(
    LocalPenaltySolveMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return R
end

function generate_model()
  model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(4,4)))
  labels = get_face_labeling(model)
  merge!(
    labels,
    Geometry.face_labeling_from_cell_tags(get_grid_topology(model),[1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3],["A","B","C"])
  )
  Geometry.add_tag_from_tags!(labels,"D",["B","C"])
  return model
end

function compute_number_of_cells(model,domain_tags,integration_tags)
  D = num_cell_dims(model)
  labels = get_face_labeling(model)

  if isnothing(domain_tags)
    domain_cell_mask = fill(true,num_cells(model))
    domain_face_mask = fill(true,num_faces(model,D-1))
  else
    domain_cell_mask = get_face_mask(labels,domain_tags,D)
    domain_face_mask = get_face_mask(labels,domain_tags,D-1)
  end
  if isnothing(integration_tags)
    integration_cell_mask = fill(true,num_cells(model))
    integration_face_mask = fill(true,num_faces(model,D-1))
  else
    integration_cell_mask = get_face_mask(labels,integration_tags,D)
    integration_face_mask = get_face_mask(labels,integration_tags,D-1)
  end
  n_domain_cells = sum(domain_cell_mask)
  n_integration_cells = sum(integration_cell_mask .& domain_cell_mask)
  n_integration_faces = sum(integration_face_mask .& domain_face_mask)

  return n_domain_cells, n_integration_cells, n_integration_faces
end

function test_operator(P,V)
  u, v = get_trial_fe_basis(V), get_fe_basis(V)
  Pu, Pv = P(u), P(v)

  uh = zero(V)
  Ruh = P(uh)
end

function test_patch_assembly(model,domain_tags,integration_tags;order=1)
  D = num_cell_dims(model)
  Ω = Triangulation(ReferenceFE{D},model;tags=domain_tags)
  Γ = Triangulation(ReferenceFE{D-1},model;tags=domain_tags)

  n_domain_cells, n_integration_cells, n_integration_faces = compute_number_of_cells(model,domain_tags,integration_tags)
  
  ptopo = Geometry.PatchTopology(model;tags=domain_tags)
  n_patches = Geometry.num_patches(ptopo)
  @test n_patches == n_domain_cells
  
  Ωp = Geometry.PatchTriangulation(model,ptopo;tags=integration_tags)
  Γp = Geometry.PatchBoundaryTriangulation(model,ptopo;tags=integration_tags)
  @test Geometry.num_cells(Ωp) == n_integration_cells
  @test length(unique(get_glue(Γp,Val(1)).tface_to_mface)) == n_integration_faces
  
  reffe = ReferenceFE(lagrangian, Float64, order; space=:P)
  VΩ = FESpace(Ω,reffe; conformity=:L2)
  VΓ = FESpace(Γ,reffe; conformity=:L2)
  X  = MultiFieldFESpace([VΩ,VΓ])
  Xb = MultiFieldFESpace([VΩ,VΓ],style=MultiField.BlockMultiFieldStyle())
  
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
  mats1 = assemble_matrix(a1,passem1,VΩ,VΩ)
  vecs1 = assemble_vector(l1,passem1,VΩ)
  matvecs1 = assemble_matrix_and_vector(a1,l1,passem1,VΩ,VΩ)
  @test length(collect(mats1)) == length(collect(vecs1)) == length(collect(matvecs1)) == n_patches
  
  passem2 = FESpaces.PatchAssembler(ptopo,VΓ,VΓ)
  a2(uΓ,vΓ) = mass(uΓ,vΓ,Γp,dΓp)
  l2(vΓ) = ∫(f⋅vΓ)dΓp
  mats2 = assemble_matrix(a2,passem2,VΓ,VΓ)
  vecs2 = assemble_vector(l2,passem2,VΓ)
  matvecs2 = assemble_matrix_and_vector(a2,l2,passem2,VΓ,VΓ,zero(VΓ))
  @test length(collect(mats2)) == length(collect(vecs2)) == length(collect(matvecs2)) == n_patches
  
  passem3 = FESpaces.PatchAssembler(ptopo,X,X)
  a3((uΩ,uΓ),(vΩ,vΓ)) = laplacian(uΩ,vΩ,dΩp) + mass(uΓ,vΩ,Γp,dΓp) + mass(uΩ,vΓ,Γp,dΓp) + mass(uΓ,vΓ,Γp,dΓp)
  l3((vΩ,vΓ)) = ∫(f⋅vΩ)dΩp
  mats3 = assemble_matrix(a3,passem3,X,X)
  vecs3 = assemble_vector(l3,passem3,X)
  matvecs3 = assemble_matrix_and_vector(a3,l3,passem3,X,X)
  @test length(collect(mats3)) == length(collect(vecs3)) == length(collect(matvecs3)) == n_patches

  passem4 = FESpaces.PatchAssembler(ptopo,Xb,Xb)
  a4((uΩ,uΓ),(vΩ,vΓ)) = laplacian(uΩ,vΩ,dΩp) + mass(uΓ,vΩ,Γp,dΓp) + mass(uΩ,vΓ,Γp,dΓp) + mass(uΓ,vΓ,Γp,dΓp)
  l4((vΩ,vΓ)) = ∫(f⋅vΩ)dΩp
  mats4 = assemble_matrix(a4,passem4,Xb,Xb)
  vecs4 = assemble_vector(l4,passem4,Xb)
  matvecs4 = assemble_matrix_and_vector(a4,l4,passem4,Xb,Xb)
  @test length(collect(mats4)) == length(collect(vecs4)) == length(collect(matvecs4)) == n_patches

  # Local Operators
  PΩ = projection_operator(VΩ,Ωp,dΩp)
  test_operator(PΩ,VΩ)

  PΓ = projection_operator(VΓ,Γp,dΓp)
  test_operator(PΓ,VΓ)
  test_operator(PΓ,VΩ)

  if n_domain_cells === n_integration_cells
    R  = reconstruction_operator(ptopo,Xb,Ωp,Γp,dΩp,dΓp,order)
    test_operator(R,Xb)

    function a_R(u,v)
      RuΩ, RuΓ = R(u)
      RvΩ, RvΓ = R(v)
      return laplacian(RuΩ,RvΩ,dΩp) + laplacian(RuΓ,RvΩ,dΩp) + laplacian(RuΩ,RvΓ,dΩp) + mass(RuΓ,RvΓ,Γp,dΓp)
    end

    Xp = FESpaces.PatchFESpace(X,ptopo)
    passem_R = FESpaces.PatchAssembler(ptopo,X,X)
    u, v = get_trial_fe_basis(X), get_fe_basis(X)
    cell_mats = FESpaces.collect_patch_cell_matrix(passem_R,Xp,Xp,a_R(u,v))
    cell_R = assemble_matrix(passem_R,cell_mats)
    @test length(collect(cell_R)) == n_patches
  end

  # Star-patch assembly 

  passem4 = FESpaces.PatchAssembler(ptopo,X,X;assembly=:star)
  a4((uΩ,uΓ),(vΩ,vΓ)) = laplacian(uΩ,vΩ,dΩp) + mass(uΓ,vΩ,Γp,dΓp) + mass(uΩ,vΓ,Γp,dΓp) + mass(uΓ,vΓ,Γp,dΓp)
  l4((vΩ,vΓ)) = ∫(f⋅vΩ)dΩp
  mats4 = assemble_matrix(a4,passem4,X,X)

  return true
end

model = generate_model()

# Case 0: 
#  - PatchTopology spans all cells
#  - Integration on all patches
test_patch_assembly(model,nothing,nothing)

# Case 1:
#  - PatchTopology spans all cells
#  - Integration on a subset of patches
test_patch_assembly(model,nothing,"D")

# Case 2:
#  - PatchTopology spans a subset of cells
#  - Integration on all patches
test_patch_assembly(model,"D",nothing)

# Case 3:
#  - PatchTopology spans a subset of cells
#  - Integration on a subset of patches
test_patch_assembly(model,"D","B")

end