using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs

##############################################################
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
##############################################################
# Toi discuss: Optimization of reconstruction operator
# The reconstruction operator matrices depend on the geometry of the element.
# If we only have uniform quads or uniform triangles then we don not need to store all the matrices.
# To check: the case of voronoi.
##############################################################

##############################################################
# To discuss: Storing unecessay zeros vs the following:
# bulk_proj_lhs(u,vT) = ∫(u*vT)dΩp  
# bulk_proj_rhs(w,vT) = ∫(w*vT)dΩp

# assem_bulk_proj_lhs = FESpaces.PatchAssembler(ptopo,X[1],Y[1])
# assem_bulk_proj_rhs = FESpaces.PatchAssembler(ptopo,Vbulk,Y[1])
# bulk_proj_lhs_mats = assemble_matrix(bulk_proj_lhs,assem_bulk_proj_lhs,X[1],Y[1])
# bulk_proj_rhs_mats = assemble_matrix(bulk_proj_rhs,assem_bulk_proj_rhs,Vbulk,Y[1])

# Π(v) = change_domain(v,Γp,DomainStyle(v))
# skel_proj_lhs(u,vF) = ∫(u*vF)dΓp 
# skel_proj_rhs(w,vF) = ∫(Π(w)*vF)dΓp

# assem_skel_proj_lhs = FESpaces.PatchAssembler(ptopo,X[2],Y[2])
# assem_skel_proj_rhs = FESpaces.PatchAssembler(ptopo,Vbulk,Y[2])
# skel_proj_lhs_mats = assemble_matrix(skel_proj_lhs,assem_skel_proj_lhs,X[2],Y[2])
# skel_proj_rhs_mats = assemble_matrix(skel_proj_rhs,assem_skel_proj_rhs,Vbulk,Y[2])
#############################################################

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
  rec_rhs((uT,uF),(v,μ)) =  ∫( (∇(uT)⋅∇(v)) + (uT*μ) )dΩp + ∫( (uF - Π(uT))*(Πn(v)) )dΓp    # (uT,uF) ∈ X, v ∈ Vr_test

  assem_rec_lhs = FESpaces.PatchAssembler(ptopo,Xr,Yr)
  assem_rec_rhs = FESpaces.PatchAssembler(ptopo,X,Yr)
  rec_lhs_mats = assemble_matrix(rec_lhs,assem_rec_lhs,Xr,Yr)
  rec_rhs_mats = assemble_matrix(rec_rhs,assem_rec_rhs,X,Yr)
  rec_op_mats = collect(lazy_map(FESpaces.HHO_ReconstructionOperatorMap(),rec_lhs_mats,rec_rhs_mats)) 

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

# projection operator
function projection_operator(X, Y, ptopo, order, Ω, dΩp, dΓp)
  reffe = ReferenceFE(lagrangian, Float64, order+1; space=:P)
  Vbulk_test= TestFESpace(Ω, reffe; conformity=:L2)
  Vbulk     = TrialFESpace(Vbulk_test)
  mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
  Xr = MultiFieldFESpace([Vbulk, X[2]];style=mfs)

  Γp = dΓp.quad.trian
  Π(v) = change_domain(v,Γp,DomainStyle(v))
  proj_lhs( (uT,uF),(vT,vF) ) = ∫(uT*vT)dΩp + ∫(uF*vF)dΓp 
  proj_rhs( (wT,wF),(vT,vF) ) = ∫(wT*vT)dΩp + ∫(Π(wT)*vF)dΓp

  assem_proj_lhs = FESpaces.PatchAssembler(ptopo,X,Y)
  assem_proj_rhs = FESpaces.PatchAssembler(ptopo,Xr,Y)

  proj_lhs_mats = assemble_matrix(proj_lhs,assem_proj_lhs,X,Y)
  proj_rhs_mats = assemble_matrix(proj_rhs,assem_proj_rhs,Xr,Y)
  proj_op_mats  = collect(lazy_map(FESpaces.HHO_L2ProjectionOperatorMap(), proj_lhs_mats, proj_rhs_mats)) 

  bulk_fields = CellData.get_data(get_fe_basis(Y[1]))
  bulk_coeffs = map(first,proj_op_mats)
  bulk_cell_basis = lazy_map(linear_combination,bulk_coeffs,bulk_fields)

  ################################################################################################
  ############ To discuss: we need a cell-wise facet linearcombination field here ################
  ################################################################################################
  skel_fields = CellData.get_data(get_fe_basis(Y[2])) 
  cell_wise_facet_dofs = Gridap.FESpaces.get_patch_dofs(Y[2],ptopo)
  cell_wise_skel_fields = map(cell_wise_facet_dofs) do dofs
    values = vcat(skel_fields...)
    cell_wise_skel_fields = Gridap.Arrays.CompressedArray(values, dofs)
    return cell_wise_skel_fields
  end
  skel_coeffs = map(last,proj_op_mats)
  skel_cell_basis = lazy_map(linear_combination,skel_coeffs,cell_wise_skel_fields)

  Proj_vΩ = FESpaces.SingleFieldFEBasis(bulk_cell_basis, Ω, FESpaces.TestBasis(), FESpaces.ReferenceDomain())
  Proj_uΩ = FESpaces.SingleFieldFEBasis(lazy_map(transpose,bulk_cell_basis), Ω, FESpaces.TrialBasis(), FESpaces.ReferenceDomain())
  ################################################### To discuss: Triangulation Γ ####################################################
  Proj_vΓ = FESpaces.SingleFieldFEBasis(skel_cell_basis, Γ, FESpaces.TestBasis(), FESpaces.ReferenceDomain()) # To discuss
  Proj_uΓ = FESpaces.SingleFieldFEBasis(lazy_map(transpose,skel_cell_basis), Γ, FESpaces.TrialBasis(), FESpaces.ReferenceDomain()) # To discuss

  return Proj_vΩ, Proj_uΩ, Proj_vΓ, Proj_uΓ
end

u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
q(x) = -∇(u)(x)
f(x) = (∇ ⋅ q)(x)

nc = (16,16)
model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),nc))
D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D}, model)
Γ = Triangulation(ReferenceFE{D-1}, model)

ptopo = Geometry.PatchTopology(model)
Ωp = Geometry.PatchTriangulation(model,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)

# Reference FEs
order = 4
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
Ybd = MultiFieldFESpace([V_test, Mbd_test];style=mfs) # For assembly and imposing boundary conditions
Xbd = MultiFieldFESpace([V, Mbd];style=mfs) # For assembly and imposing boundary conditions
Y   = MultiFieldFESpace([V_test, M_test];style=mfs) # For computing local opearotrs
X   = MultiFieldFESpace([V, M];style=mfs) # For computing local opearotrs

degree = order+2 
dΩp = Measure(Ωp,degree)
dΓp = Measure(Γp,degree)

R_vΩ, R_uΩ, R_vΓ, R_uΓ  = reconstruction_operator(X, ptopo, Ω, dΩp, dΓp, order)
P_vΩ, P_uΩ, P_vΓ, P_uΓ  = projection_operator(X, Y, ptopo, order, Ω, dΩp, dΓp)



