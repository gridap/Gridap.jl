using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs


# To discuss:
# 1. get_patch_dofs for zeromean fe spaces  
# # function get_patch_dofs(space::FESpace,ptopo::PatchTopology)
#   trian = get_triangulation(space)

#   Df = num_cell_dims(trian)
#   face_to_tface = get_glue(trian,Val(Df)).mface_to_tface
#   @notimplementedif isnothing(face_to_tface)
#   tface_to_dofs = get_cell_dof_ids(space)
#   patch_to_faces = Gridap.FESpaces.get_patch_faces(ptopo,Df)

#   using DataStructures: SortedSet
#   using Gridap.Arrays: Table


## To discuss: The bug for zeromean fe spaces occurs here ##
#   patch_to_dofs = map(patch_to_faces) do pfaces
#     tfaces = filter(x->x>0,face_to_tface[pfaces])
#     dofs = SortedSet{Int}()
#     for tface in tfaces
#       tface_dofs = filter(x->x>0,view(tface_to_dofs,tface))
#       !isempty(tface_dofs) && push!(dofs,tface_dofs...)
#     end
#     collect(dofs)
#   end |> Table
##############################################################
# #   return patch_to_dofs
# # end

# Toi discuss: Optimization of reconstruction operator
# The reconstruction operator matrices depend on the geometry of the element.
# If we only have uniform quads or uniform triangles then we don not need to store all the matrices.
# To check: the case of voronoi.

function reconstruction_operator(X, ptopo, Ω, dΩp, dΓp, order)

  Ωp = dΩp.quad.trian
  Γp = dΓp.quad.trian

  reffeVrec = ReferenceFE(lagrangian, Float64, order+1; space=:P)
  reffeC    = ReferenceFE(lagrangian, Float64, 0; space=:P)
  
  # Vr_test  = TestFESpace(Ω, reffeVr; conformity=:L2, constraint=:zeromean) 
  Vrec_test = TestFESpace(Ω, reffeVrec; conformity=:L2)
  C_test  = TestFESpace(Ω, reffeC; conformity=:L2)   
  Vrec = TrialFESpace(Vrec_test) 
  C  = TrialFESpace(C_test)

  mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
  Yr = MultiFieldFESpace([Vrec_test, C_test];style=mfs)
  Xr = MultiFieldFESpace([Vrec, C];style=mfs)

  nrel = get_normal_vector(Γp)
  Πn(v) = ∇(v)⋅nrel
  Π(v) = change_domain(v,Γp,DomainStyle(v))

  rec_lhs((u,λ),(v,μ))   = ∫( (∇(u)⋅∇(v)) + (μ*u) + (λ*v) )dΩp   # u ∈ Vr_trial, v ∈ Vr_test
  rec_rhs((uT,uF),(v,μ)) =  ∫( (∇(uT)⋅∇(v)) + (uT*μ) )dΩp + ∫( (Π(uF) - uT)*(Πn(v)) )dΓp    # (uT,uF) ∈ X, v ∈ Vr_test

  assem_rec_lhs = FESpaces.PatchAssembler(ptopo,Xr,Yr)
  assem_rec_rhs = FESpaces.PatchAssembler(ptopo,X,Yr)
  rec_lhs_mats = assemble_matrix(rec_lhs,assem_rec_lhs,Xr,Yr)
  rec_rhs_mats = assemble_matrix(rec_rhs,assem_rec_rhs,X,Yr)
  rec_op_mats = collect(lazy_map(FESpaces.ReconstructionOperatorMap(),rec_lhs_mats,rec_rhs_mats)) 

  fields = CellData.get_data(get_fe_basis(Vrec_test))
  coeffsΩ = map(first,rec_op_mats)
  cell_basis_Ω = lazy_map(linear_combination,coeffsΩ,fields)
  coeffsΓ = map(last,rec_op_mats)
  cell_basis_Γ = lazy_map(linear_combination,coeffsΓ,fields)

  rec_vΩ = FESpaces.SingleFieldFEBasis(cell_basis_Ω, Ω, FESpaces.TestBasis(), FESpaces.ReferenceDomain())
  rec_uΩ = FESpaces.SingleFieldFEBasis(lazy_map(transpose,cell_basis_Ω), Ω, FESpaces.TrialBasis(), FESpaces.ReferenceDomain())
  rec_vΓ = FESpaces.SingleFieldFEBasis(cell_basis_Γ, Ω, FESpaces.TestBasis(), FESpaces.ReferenceDomain())
  rec_uΓ = FESpaces.SingleFieldFEBasis(lazy_map(transpose,cell_basis_Γ), Ω, FESpaces.TrialBasis(), FESpaces.ReferenceDomain())

  return rec_vΩ, rec_uΩ, rec_vΓ, rec_uΓ 
end

function statically_condensed_assembly(ptopo,X,Y,M,M_test,a,l)
  # Lazily assemble the global patch-systems and 
  # perform the static condensation
  assem = FESpaces.PatchAssembler(ptopo,X,Y)
  full_matvecs = assemble_matrix_and_vector(a,l,assem,X,Y)
  sc_matvecs = lazy_map(FESpaces.StaticCondensationMap(),full_matvecs)

  # Regular assembly of the statically-assembled systems
  patch_rows = assem.strategy.array.array[2,2].patch_rows
  patch_cols = assem.strategy.array.array[2,2].patch_cols
  matvecdata = ([sc_matvecs,],[patch_rows,],[patch_cols,])
  matdata = ([],[],[]) # dummy matdata
  vecdata = ([],[],[]) # dummy vecdata
  data = (matvecdata, matdata, vecdata)
  A, b = assemble_matrix_and_vector(SparseMatrixAssembler(M,M_test),data)
  return A, b
end

function backward_static_condensation(ptopo,X,Y,a,l,xb)
  assem = FESpaces.PatchAssembler(ptopo,X,Y)
  patch_rows_bb = assem.strategy.array.array[2,2].patch_rows
  patch_cols_bb = assem.strategy.array.array[2,2].patch_cols

  # To discuss: can patch_rows and patch_cols be different from each other for a patch?
  if patch_rows_bb != patch_cols_bb
    @notimplemented
  else
    patchwise_xb = map(patch_rows_bb) do patch_row
      patchwise_xb = xb[patch_row]
      return patchwise_xb
    end
  end

  full_matvecs = assemble_matrix_and_vector(a,l,assem,X,Y)
  patchwise_xi_vec = lazy_map(Gridap.FESpaces.BackwardStaticCondensationMap(),full_matvecs,patchwise_xb)

  patch_rows_ii = assem.strategy.array.array[1,1].patch_rows
  patch_cols_ii = assem.strategy.array.array[1,1].patch_cols
  patchwise_xi_vecdata = ([patchwise_xi_vec,],[patch_rows_ii,],[patch_cols_ii,])

  Xi = MultiFieldFESpace([X[1],X[2]])
  Yi = MultiFieldFESpace([Y[1],Y[2]])
  assem = SparseMatrixAssembler(Xi,Yi)
  xi = assemble_vector(assem, patchwise_xi_vecdata)
  wh = FEFunction(Yi,xi)
  return wh
end

u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
q(x) = -∇(u)(x)
f(x) = (∇ ⋅ q)(x)

nc = (2,2)
model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),nc))
D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D}, model)
Γ = Triangulation(ReferenceFE{D-1}, model)

ptopo = Geometry.PatchTopology(model)
Ωp = Geometry.PatchTriangulation(model,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)

# Reference FEs
order = 0
reffeV = ReferenceFE(lagrangian, Float64, order; space=:P)
reffeM = ReferenceFE(lagrangian, Float64, order; space=:P)

# HHO test FE Spaces
V_test = TestFESpace(Ω, reffeV; conformity=:L2)
Mbd_test = TestFESpace(Γ, reffeM; conformity=:L2, dirichlet_tags="boundary")  # For assembly and imposing boundary conditions
M_test = TestFESpace(Γ, reffeM; conformity=:L2)   # For computing local opearotrs 

# HDG trial FE Spaces
V   = TrialFESpace(V_test)
M   = TrialFESpace(M_test) 
Mbd = TrialFESpace(Mbd_test, u)

mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
Ybd = MultiFieldFESpace([V_test, Mbd_test];style=mfs)
Xbd = MultiFieldFESpace([V, Mbd];style=mfs)
Y   = MultiFieldFESpace([V_test, M_test];style=mfs)
X   = MultiFieldFESpace([V, M];style=mfs)

degree = order+2 
dΩp = Measure(Ωp,degree)
dΓp = Measure(Γp,degree)

R_vΩ, R_uΩ, R_vΓ, R_uΓ  = reconstruction_operator(X, ptopo, Ω, dΩp, dΓp, order)






