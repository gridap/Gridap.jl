using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs, Gridap.Arrays


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

  nfields = length(X.spaces)
  RvΩ = MultiField.MultiFieldFEBasisComponent(rec_vΩ,1,nfields)
  RuΩ = MultiField.MultiFieldFEBasisComponent(rec_uΩ,1,nfields)
  RvΓ = MultiField.MultiFieldFEBasisComponent(rec_vΓ,2,nfields)
  RuΓ = MultiField.MultiFieldFEBasisComponent(rec_uΓ,2,nfields)

  return RvΩ, RuΩ, RvΓ, RuΓ 
end

function projection_operator(u, order, Ωs, Ωp)
  reffe = ReferenceFE(lagrangian, Float64, order; space=:P)
  P = FESpace(Ωs, reffe; conformity=:L2)

  dΩp = Measure(Ωp,2*(order+1))
  mass(u,v) = ∫(u⋅v)*dΩp
  
  uP = get_trial_fe_basis(P)
  vP = change_domain(get_fe_basis(P),Ωp,FESpaces.ReferenceDomain())
  lhs = get_contribution(mass(uP,vP),Ωp)
  rhs = get_contribution(mass(u,vP),Ωp)
  fields = CellData.get_data(vP)
  coeffs = lazy_map(FESpaces.LocalSolveMap(),lhs,rhs)
  cell_basis = lazy_map(linear_combination,coeffs,fields)
  proj_v = FESpaces.SingleFieldFEBasis(cell_basis, Ωp, FESpaces.TestBasis(), FESpaces.ReferenceDomain())
  Proj_u = FESpaces.SingleFieldFEBasis(lazy_map(transpose, cell_basis), Ωp, FESpaces.TrialBasis(), FESpaces.ReferenceDomain())
  return proj_v, Proj_u
end


function mixed_order_projection(X, Y, order, Γ, Γp)
  Π(v) = change_domain(v,Γp,DomainStyle(v))
  uT = get_trial_fe_basis(X[1]);
  PΓ_vT, PΓ_uT = projection_operator(Π(uT), order, Γ, Γp)

  nfields = length(X.spaces)
  PΓ_vT =  MultiField.MultiFieldFEBasisComponent(PΓ_vT, 1, nfields)
  PΓ_uT =  MultiField.MultiFieldFEBasisComponent(PΓ_uT, 1, nfields)

  vT, vF = get_fe_basis(Y);
  uT, uF = get_trial_fe_basis(X);

  δTFv = PΓ_vT - Π(vF)
  δTFu = PΓ_uT - Π(uF)

  return δTFv, δTFu
end

function statically_condensed_assembly(ptopo,X,Y,Mbd,Mbd_test,a_consistency,a_stabilization,load)

  # patchwise assembly 
  # If I try to assemble with Xb and Yb here. Then the boundary information is not passed to the consistency matrix.
  assemp = FESpaces.PatchAssembler(ptopo,X,Y)
  ############ ac_cellwise_mats is Probably incorrect, since get_cell_dof_ids(test, trian) is incorrect ################## 
  ac_cellwise_mats, _, _ = collect_cell_matrix(X,Y,a_consistency)
  ########################################################################################################################
  as_patchwise_mats = assemble_matrix(assemp, FESpaces.collect_patch_cell_matrix(assemp,X,Y,a_stabilization))
  l_patchwise_vecs  = assemble_vector(assemp, FESpaces.collect_patch_cell_vector(assemp,Y,load))

  # Not able to collect collect lhs mats (consistency and stabilization) together. So proceeding as follows.
  ac_cellwise_mats = first(ac_cellwise_mats)
  sc_matvecs = map(FESpaces.HHOStaticCondensationMap(), ac_cellwise_mats, as_patchwise_mats,l_patchwise_vecs)

  assem = FESpaces.PatchAssembler(ptopo,Xbd,Ybd)
  patch_rows = assem.strategy.array.array[2,2].patch_rows
  patch_cols = assem.strategy.array.array[2,2].patch_cols
  matvecdata = ([sc_matvecs,],[patch_rows,],[patch_cols,])
  matdata = ([],[],[]) # dummy matdata
  vecdata = ([],[],[]) # dummy vecdata
  data = (matvecdata, matdata, vecdata)
  Asc, bsc = assemble_matrix_and_vector(SparseMatrixAssembler(Mbd,Mbd_test),data) # assemble and apply boundary conditions
  return Asc, bsc
end


##############################################################
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
reffeV = ReferenceFE(lagrangian, Float64, order+1; space=:P)
reffeM = ReferenceFE(lagrangian, Float64, order; space=:P)

# HHO test FE Spaces
V_test = TestFESpace(Ω, reffeV; conformity=:L2)
M_test = TestFESpace(Γ, reffeM; conformity=:L2)   # For computing local opearotrs 
Mbd_test = TestFESpace(Γ, reffeM; conformity=:L2, dirichlet_tags="boundary")  # For assembly and imposing boundary conditions

V   = TrialFESpace(V_test)
M   = TrialFESpace(M_test) 
Mbd = TrialFESpace(Mbd_test, u)

mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
Ybd = MultiFieldFESpace([V_test, Mbd_test];style=mfs) # For assembly and imposing boundary conditions
Xbd = MultiFieldFESpace([V, Mbd];style=mfs) # For assembly and imposing boundary conditions
Y   = MultiFieldFESpace([V_test, M_test];style=mfs) # For computing local opearotrs
X   = MultiFieldFESpace([V, M];style=mfs) # For computing local opearotrs

degree = 2*(order+1) 
dΩp = Measure(Ωp,degree)
dΓp = Measure(Γp,degree)

RvΩ, RuΩ, RvΓ, RuΓ  = reconstruction_operator(X, ptopo, Ω, dΩp, dΓp, order)

hT=CellField(get_array(∫(1)dΩp),Ωp)
hTinv = 1/hT
δTFv, δTFu = mixed_order_projection(X, Y, order, Γ, Γp)

a_consistency   = ∫( (RuΩ*RvΩ)+(RuΩ*RvΓ)+(RuΓ*RvΩ)+(RuΓ*RvΓ) )dΩp 
a_stabilization = ∫( hTinv*(δTFv*δTFu) )dΓp

vT, vF = get_fe_basis(Y);
load = ∫( f*vT )*dΩp

Asc, bsc = statically_condensed_assembly(ptopo,X,Y,Mbd,Mbd_test,a_consistency,a_stabilization,load)

solver = LUSolver()
ns = numerical_setup(symbolic_setup(solver,Asc),Asc)
xb = zeros(size(bsc))
solve!(xb,ns,bsc)
sh = FEFunction(Mbd_test,xb)

