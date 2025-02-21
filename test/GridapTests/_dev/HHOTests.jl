using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs


##############################################################
# Toi discuss: Optimization of reconstruction operator
# The reconstruction operator matrices depend on the geometry of the element.
# If we only have uniform quads or uniform triangles then we don not need to store all the matrices.
# To check: the case of voronoi.
##############################################################


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

function projection_operator(u, order, Ωs, Ωp)
  reffe = ReferenceFE(lagrangian, Float64, order+1; space=:P)
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

function consistency(rec_vΩ, rec_uΩ, rec_vΓ, rec_uΓ, X)

  nfields = length(X.spaces)
  RvΩ = MultiField.MultiFieldFEBasisComponent(rec_vΩ,1,nfields)
  RuΩ = MultiField.MultiFieldFEBasisComponent(rec_uΩ,1,nfields)
  RvΓ = MultiField.MultiFieldFEBasisComponent(rec_vΓ,2,nfields)
  RuΓ = MultiField.MultiFieldFEBasisComponent(rec_uΓ,2,nfields)

    return ∫( (RuΩ*RvΩ)+(RuΩ*RvΓ)+(RuΓ*RvΩ)+(RuΓ*RvΓ) )Measure(Ω,4)
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
reffeV = ReferenceFE(lagrangian, Float64, order; space=:P)
reffeM = ReferenceFE(lagrangian, Float64, order; space=:P)

# HHO test FE Spaces
V_test = TestFESpace(Ω, reffeV; conformity=:L2)
Mbd_test = TestFESpace(Γ, reffeM; conformity=:L2, dirichlet_tags="boundary")  # For assembly and imposing boundary conditions
M_test = TestFESpace(Γ, reffeM; conformity=:L2)   # For computing local opearotrs 

V   = TrialFESpace(V_test)
M   = TrialFESpace(M_test) 
Mbd = TrialFESpace(Mbd_test, u)

mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
Ybd = MultiFieldFESpace([V_test, Mbd_test];style=mfs) # For assembly and imposing boundary conditions
Xbd = MultiFieldFESpace([V, Mbd];style=mfs) # For assembly and imposing boundary conditions
Y   = MultiFieldFESpace([V_test, M_test];style=mfs) # For computing local opearotrs
X   = MultiFieldFESpace([V, M];style=mfs) # For computing local opearotrs

Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
dΩ = Measure(Ω,2*(order+1))
mass(u,v,dΩ) = ∫(u⋅v)dΩ
massΩ(u,v) = mass(u,v,dΩ)
massΓ(u,v) = mass(u,Π(v,Γp),dΓp)
PΩ = FESpaces.LocalOperator(
  FESpaces.LocalSolveMap(), V, massΩ, massΩ
)
PΓ = FESpaces.LocalOperator(
  FESpaces.LocalSolveMap(), M, massΓ, massΓ; trian_out = Γp
)

R_reffe = ReferenceFE(lagrangian, Float64, order+1; space=:P)
R_space = FESpace(Ω, R_reffe; conformity=:L2)
vR = get_fe_basis(R_space)
PvΩ = evaluate(PΩ,vR)
PvΓ = evaluate(PΓ,vR)

degree = order+2 
dΩp = Measure(Ωp,degree)
dΓp = Measure(Γp,degree)

RvΩ, RuΩ, RvΓ, RuΓ  = reconstruction_operator(X, ptopo, Ω, dΩp, dΓp, order)

# V = FESpace(Ω, ReferenceFE(lagrangian,Float64,order+1;space=:P); conformity=:L2)
# uh = get_fe_basis(V)
# PvΩ  = projection_operator(uh, order, Ω, Ωp)
# PvΓ  = projection_operator(uh, order, Γ, Γp)

a_consistency = consistency(RvΩ, RuΩ, RvΓ, RuΓ, X) 

assem = FESpaces.PatchAssembler(ptopo,X,Y)
cell_mats = collect_cell_matrix(X,Y,a_consistency)
# full_mats = assemble_matrix(SparseMatrixAssembler(X,Y),cell_mats)
patch_rows_ii = assem.strategy.array.array[1,1].patch_rows
# patch_col_ii = assem.strategy.array.array[1,1].patch_cols
patch_rows_ib = assem.strategy.array.array[1,2].patch_rows
# patch_col_ib = assem.strategy.array.array[1,2].patch_cols
patch_rows_bi = assem.strategy.array.array[2,1].patch_rows
patch_rows_bb = assem.strategy.array.array[2,2].patch_rows



# matvecdata = ([cell_mats,],[patch_rows,],[patch_cols,])
# matdata = ([],[],[]) # dummy matdata
# # vecdata = ([],[],[]) # dummy vecdata
# data = (matvecdata, matdata, vecdata)
# A, b = assemble_matrix_and_vector(SparseMatrixAssembler(M,M_test),data)


# full_matvecs = assemble_matrix(a_consistency,assem,Xbd,Ybd)
# sc_matvecs = lazy_map(FESpaces.StaticCondensationMap(),full_matvecs)

# # Regular assembly of the statically-assembled systems
# patch_rows = assem.strategy.array.array[2,2].patch_rows
# patch_cols = assem.strategy.array.array[2,2].patch_cols
# matvecdata = ([sc_matvecs,],[patch_rows,],[patch_cols,])
# matdata = ([],[],[]) # dummy matdata
# vecdata = ([],[],[]) # dummy vecdata
# data = (matvecdata, matdata, vecdata)
# A, b = assemble_matrix_and_vector(SparseMatrixAssembler(M,M_test),data)

hF = CellField(get_array(∫(1)dΓp),Γp)
hFinv = 1/hF

PvΩ_RΩ, PuΩ_RΩ = projection_operator(RuΩ, order, Ω, Ωp)
PvΩ_RΓ, PuΩ_RΓ = projection_operator(RuΓ, order, Ω, Ωp)
PvΓ_RΩ, PuΓ_RΩ = projection_operator(RuΩ, order, Γ, Γp)
PvΓ_RΓ, PuΓ_RΓ = projection_operator(RuΓ, order, Γ, Γp)

nfields = length(X.spaces)

mf_PuΩ_RΩ = MultiField.MultiFieldFEBasisComponent(PuΩ_RΩ,1,nfields)
mf_PuΩ_RΓ = MultiField.MultiFieldFEBasisComponent(PuΩ_RΓ,2,nfields)
mf_PuΓ_RΩ = MultiField.MultiFieldFEBasisComponent(PuΓ_RΩ,1,nfields)
mf_PuΓ_RΓ = MultiField.MultiFieldFEBasisComponent(PuΓ_RΓ,2,nfields)

mf_PvΩ_RΩ = MultiField.MultiFieldFEBasisComponent(PvΩ_RΩ,1,nfields)
mf_PvΩ_RΓ = MultiField.MultiFieldFEBasisComponent(PvΩ_RΓ,2,nfields)
mf_PvΓ_RΩ = MultiField.MultiFieldFEBasisComponent(PvΓ_RΩ,1,nfields)
mf_PvΓ_RΓ = MultiField.MultiFieldFEBasisComponent(PvΓ_RΓ,2,nfields)

Π(v) = change_domain(v,Γp,DomainStyle(v))

uT, uF = get_trial_fe_basis(X);
vT, vF = get_fe_basis(Y);

δFT = ∫(hFinv * (uF - Π(uT))*(vF - Π(vT)) )*dΓp
PΩ_R = ∫(hFinv * (mf_PuΩ_RΩ)*(mf_PvΩ_RΩ) )*dΓp + ∫(hFinv * (mf_PuΩ_RΓ)*(mf_PvΩ_RΓ) )*dΓp  # Cannot club them together. strian !=== ttrian. Why?
PΓ_R = ∫(hFinv * (mf_PuΓ_RΩ)*(mf_PvΓ_RΩ) )*dΓp + ∫(hFinv * (mf_PuΓ_RΓ)*(mf_PvΓ_RΓ) )*dΓp # Cannot club them together. strian !=== ttrian. Why?

arr = get_array(PΩ_R - PΓ_R)
 
δFT_arr = get_array(δFT)
PΩ_R_arr = get_array(PΩ_R);
PΓ_R_arr = get_array(PΓ_R);


