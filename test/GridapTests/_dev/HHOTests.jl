using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs

function consistency(rec_vΩ, rec_uΩ, rec_vΓ, rec_uΓ, X)

  nfields = length(X.spaces)
  RvΩ = MultiField.MultiFieldFEBasisComponent(rec_vΩ,1,nfields)
  RuΩ = MultiField.MultiFieldFEBasisComponent(rec_uΩ,1,nfields)
  RvΓ = MultiField.MultiFieldFEBasisComponent(rec_vΓ,2,nfields)
  RuΓ = MultiField.MultiFieldFEBasisComponent(rec_uΓ,2,nfields)

    return ∫( (RuΩ*RvΩ)+(RuΩ*RvΓ)+(RuΓ*RvΩ)+(RuΓ*RvΓ) )Measure(Ω,4)
end

function projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  P = FESpaces.LocalOperator(
    FESpaces.LocalSolveMap(), V, mass, mass
  )
  return P
end

function reconstruction_operator(ptopo,L,X,Ωp,Γp,dΩp,dΓp)
  reffe_Λ = ReferenceFE(lagrangian, Float64, 0; space=:P)
  Λ = FESpace(Ω, reffe_Λ; conformity=:L2)

  nrel = get_normal_vector(Γp)
  Πn(v) = ∇(v)⋅nrel
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((u,λ),(v,μ))   = ∫( (∇(u)⋅∇(v)) + (μ*u) + (λ*v) )dΩp
  rhs((uT,uF),(v,μ)) =  ∫( (∇(uT)⋅∇(v)) + (uT*μ) )dΩp + ∫( (uF - Π(uT,Γp))*(Πn(v)) )dΓp
  
  W = MultiFieldFESpace([L,Λ];style=mfs)
  R = FESpaces.LocalOperator(
    FESpaces.HHO_ReconstructionOperatorMap(), ptopo, W, X, lhs, rhs; space_out = L
  )
  return R
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

order = 0
qdegree = order+1

dΩ = Measure(Ω,qdegree)
dΩp = Measure(Ωp,qdegree)
dΓp = Measure(Γp,qdegree)

reffe_V = ReferenceFE(lagrangian, Float64, order; space=:P)   # Bulk space
reffe_M = ReferenceFE(lagrangian, Float64, order; space=:P)   # Skeleton space
reffe_L = ReferenceFE(lagrangian, Float64, order+1; space=:P) # Reconstruction space
V = FESpace(Ω, reffe_V; conformity=:L2)
M = FESpace(Γ, reffe_M; conformity=:L2)
L = FESpace(Ω, reffe_L; conformity=:L2)

mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
X   = MultiFieldFESpace([V, M];style=mfs)

PΩ = projection_operator(V, Ω, dΩ)
PΓ = projection_operator(M, Γp, dΓp)

R = reconstruction_operator(ptopo,L,X,Ωp,Γp,dΩp,dΓp)

Ru, Rv = evaluate(R,X)
RuΩ, RuΓ = Ru
RvΩ, RvΓ = Rv

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


Mbd = TestFESpace(Γ, reffeM; conformity=:L2, dirichlet_tags="boundary")  # For assembly and imposing boundary conditions
Mbd = TrialFESpace(Mbd_test, u)


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


