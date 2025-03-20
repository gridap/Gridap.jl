using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs

function projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  P = FESpaces.LocalOperator(
    FESpaces.LocalSolveMap(), V0, mass, mass; trian_out = Ω
  )
  return P
end

function reconstruction_operator(ptopo,L,X,Ω,Γp,dΩp,dΓp)
  reffe_Λ = ReferenceFE(lagrangian, Float64, 0; space=:P)
  Λ = FESpace(Ω, reffe_Λ; conformity=:L2)

  nrel = get_normal_vector(Γp)
  Πn(v) = ∇(v)⋅nrel
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((u,λ),(v,μ))   =  ∫( (∇(u)⋅∇(v)) + (μ*u) + (λ*v) )dΩp
  rhs((uT,uF),(v,μ)) =  ∫( (∇(uT)⋅∇(v)) + (uT*μ) )dΩp + ∫( (uF - Π(uT,Γp))*(Πn(v)) )dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  mfs = Y.multi_field_style
  W = MultiFieldFESpace([L,Λ];style=mfs)
  R = FESpaces.LocalOperator(
    FESpaces.HHO_ReconstructionOperatorMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return R
end

function projection_operator(X, Ωp, Γp, PΓ::FESpaces.LocalOperator, PΩRΩ::FESpaces.LocalOperator, PΩRΓ::FESpaces.LocalOperator)
  
  u = get_trial_fe_basis(X)
  
  RuΩ, RuΓ = R(u)
  _RuΩ = FESpaces.SingleFieldFEBasis(lazy_map(transpose,RuΩ.cell_field.args[1].args[1]),Ωp,FESpaces.TrialBasis(),ReferenceDomain())
  _RuΓ = FESpaces.SingleFieldFEBasis(lazy_map(transpose,RuΓ.cell_field.args[1].args[1]),Ωp,FESpaces.TrialBasis(),ReferenceDomain())


  PΩ_RuΩ = PΩRΩ(_RuΩ)
  PΩ_RuΓ = PΩRΓ(_RuΓ)

  wΩ = FESpaces.SingleFieldFEBasis(CellData.get_data(_RuΩ - PΩ_RuΩ),Ωp,FESpaces.TrialBasis(),ReferenceDomain())
  wΓ = FESpaces.SingleFieldFEBasis(CellData.get_data(_RuΓ - PΩ_RuΓ),Ωp,FESpaces.TrialBasis(),ReferenceDomain())

  PΓwΩ = PΓ(wΩ)
  PΓwΓ = PΓ(wΓ)  

  _PΓwΩ = FESpaces.SingleFieldFEBasis(CellData.get_data(PΓwΩ),Γp,FESpaces.TrialBasis(),ReferenceDomain())
  _PΓwΓ = FESpaces.SingleFieldFEBasis(CellData.get_data(PΓwΓ),Γp,FESpaces.TrialBasis(),ReferenceDomain())
  _PΓsΩ = FESpaces.SingleFieldFEBasis(lazy_map(transpose,CellData.get_data(PΓwΩ)),Γp,FESpaces.TestBasis(),ReferenceDomain())
  _PΓsΓ = FESpaces.SingleFieldFEBasis(lazy_map(transpose,CellData.get_data(PΓwΓ)),Γp,FESpaces.TestBasis(),ReferenceDomain())

  mf_PΓwΩ = MultiField.MultiFieldFEBasisComponent(_PΓwΩ,1,2)
  mf_PΓwΓ = MultiField.MultiFieldFEBasisComponent(_PΓwΓ,2,2)
  mf_PΓsΩ = MultiField.MultiFieldFEBasisComponent(_PΓsΩ,1,2)
  mf_PΓsΓ = MultiField.MultiFieldFEBasisComponent(_PΓsΓ,2,2)

  return mf_PΓwΩ, mf_PΓwΓ, mf_PΓsΩ, mf_PΓsΓ
end

##############################################################
uex(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
f(x) = -Δ(uex)(x)

nc = (2,2)
model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),nc))
D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D}, model)
Γ = Triangulation(ReferenceFE{D-1}, model)

ptopo = Geometry.PatchTopology(model)
Ωp = Geometry.PatchTriangulation(model,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)

order = 1
qdegree = 2*(order+1)

dΩp = Measure(Ωp,qdegree)
dΓp = Measure(Γp,qdegree)

##########################
# Mixed order variant
##########################
reffe_V = ReferenceFE(lagrangian, Float64, order; space=:P)     # Bulk space
reffe_M = ReferenceFE(lagrangian, Float64, order; space=:P)     # Skeleton space
reffe_L = ReferenceFE(lagrangian, Float64, order+1; space=:P)   # Reconstruction space
V = FESpace(Ω, reffe_V; conformity=:L2)
M = FESpace(Γ, reffe_M; conformity=:L2, dirichlet_tags="boundary")
L = FESpace(Ω, reffe_L; conformity=:L2)
N = TrialFESpace(M,uex)

mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
X  = MultiFieldFESpace([V, N];style=mfs)
Y  = MultiFieldFESpace([V, M];style=mfs)
Xp = FESpaces.PatchFESpace(X,ptopo)

R  = reconstruction_operator(ptopo,L,Y,Ωp,Γp,dΩp,dΓp)
PΓ = projection_operator(M, Γp, dΓp)


######
Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
mass(u,v,Ω,dΩ) = ∫(u⋅Π(v,Ω))dΩ
mass_Ω(u,v) = mass(u,v,Ωp,dΩp)
M0 = FESpaces.FESpaceWithoutBCs(M)
M0p = FESpaces.PatchFESpace(M0,ptopo)

PΩRΩ = FESpaces.LocalOperator(FESpaces.LocalSolveMap(), ptopo, V, V, mass_Ω, mass_Ω)
PΩRΓ = FESpaces.LocalOperator(FESpaces.LocalSolveMap(), ptopo, V, M0p, mass_Ω, mass_Ω)

# (PΓ ∘ (I - PΩ) ∘ R)(u) = (PΓ((I-PΩ)(R(u))))
PΓRuΩ, PΓRuΓ, PΓRvΩ, PΓRvΓ = projection_operator(X, Ωp, Γp, PΓ, PΩRΩ, PΩRΓ)
#####

Xp = FESpaces.PatchFESpace(X,ptopo)
Yp = FESpaces.PatchFESpace(Y,ptopo)
global_assem = SparseMatrixAssembler(X,Y)
patch_assem = FESpaces.PatchAssembler(ptopo,X,Y)

function a(u,v)
  Ru_Ω, Ru_Γ = R(u)
  Rv_Ω, Rv_Γ = R(v)
  return ∫(∇(Ru_Ω)⋅∇(Rv_Ω) + ∇(Ru_Γ)⋅∇(Rv_Ω) + ∇(Ru_Ω)⋅∇(Rv_Γ) + ∇(Ru_Γ)⋅∇(Rv_Γ))dΩp
end

hTinv =  CellField(1 ./ (sqrt(2).*sqrt.(get_array(∫(1)dΩp))),Ωp)
function s(u,v)
  function SΓa(u)
    u_Ω, u_Γ = u
    return PΓ(u_Ω) - u_Γ 
  end
  SΓa_u = SΓa(u)
  SΓa_v = SΓa(v)
  SΓb_u = (PΓRuΩ + PΓRuΓ)
  SΓb_v = (PΓRvΩ + PΓRvΓ)
  ################### Breaks ###################
  # ∫( hTinv * ( ( SΓ(u) + (PΓRuΩ + PΓRuΓ) ) * ( SΓ(v) + (PΓRvΩ + PΓRvΓ) ) ) )dΓp  #
  ##############################################
  return ∫( hTinv * ( SΓa_u * SΓa_v ) )dΓp + ∫( hTinv * (SΓa_u * SΓb_v) )dΓp + ∫( hTinv * (SΓa_v * SΓb_u) )dΓp + ∫( hTinv * (SΓb_u * SΓb_v) )dΓp
end

l((vΩ,vΓ)) = ∫(f⋅vΩ)dΩp

u, v = get_trial_fe_basis(X), get_fe_basis(Y)

function SΓa(u)
  u_Ω, u_Γ = u
  return PΓ(u_Ω) - u_Γ 
end
SΓb_u = (PΓRuΩ + PΓRuΓ)
SΓb_v = (PΓRvΩ + PΓRvΓ)

SΓa_u = SΓa(u)
SΓa_v = SΓa(v)

x = get_cell_points(Γp)

PΓRuΩ(x)[1]
PΓRvΩ(x)[1]
PΓRuΓ(x)[1]
PΓRvΓ(x)[1]

SΓb_u(x)[1]
SΓb_v(x)[1]

SΓa_u(x)[1]

mass_Γ(u,v) = mass(u,v,Γp,dΓp)

function weakform()
  ll = FESpaces.collect_cell_matrix_and_vector(Xp,Xp,a(u,v),DomainContribution(),zero(Xp))
  aa = FESpaces.collect_cell_matrix_and_vector(X,X,mass_Γ(SΓa_u,SΓa_v),l(v),zero(X))
  ab = FESpaces.collect_cell_matrix_and_vector(X,Xp,mass_Γ(SΓa_u,SΓb_v),DomainContribution(),zero(X))
  ba = FESpaces.collect_cell_matrix_and_vector(Xp,X,mass_Γ(SΓb_u,SΓa_v),DomainContribution(),zero(Xp))
  bb = FESpaces.collect_cell_matrix_and_vector(Xp,Xp,mass_Γ(SΓb_u,SΓb_v),DomainContribution(),zero(Xp))
  data = FESpaces.merge_assembly_matvec_data(ll,aa,ab,ba,bb)
  assemble_matrix_and_vector(global_assem,data)
end

A, b = weakform()
x = A \ b

uΩ, uΓ = FEFunction(X,x)

eh = uΩ - uex
sum(∫(eh⋅eh)dΩp)


∫(hTinv * (SΓb_u * SΓb_v))dΓp

S1 = ∫( hTinv * (SΓa(u) * SΓa(v)) )dΓp
S2 = ∫( hTinv * (SΓa(u) * SΓb_v) )dΓp
S3 = ∫( hTinv * (SΓa(v) * SΓb_u) )dΓp
S4 = ∫( hTinv * (SΓb_u * SΓb_v) )dΓp

data1 = FESpaces.collect_cell_matrix_and_vector(X,Y,s(u,v),l(v),zero(X))



function weakform()
  u, v = get_trial_fe_basis(X), get_fe_basis(Y)
  data1 = FESpaces.collect_cell_matrix_and_vector(X,Y,s(u,v),l(v),zero(X))
  data2 = FESpaces.collect_cell_matrix_and_vector(Xp,Xp,a(u,v),DomainContribution(),zero(Xp))
  data = FESpaces.merge_assembly_matvec_data(data1,data2)
  assemble_matrix_and_vector(global_assem,data)
end

function patch_weakform()
  u, v = get_trial_fe_basis(X), get_fe_basis(Y)
  data1 = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,X,Y,s(u,v),l(v),zero(X))
  data2 = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,Xp,Xp,a(u,v),DomainContribution(),zero(Xp))
  data = FESpaces.merge_assembly_matvec_data(data1,data2)
  return assemble_matrix_and_vector(patch_assem,data)
end

# Monolithic solve
A, b = weakform()
x = A \ b

# Static condensation
op = MultiField.StaticCondensationOperator(X,V,N,patch_assem,patch_weakform())

ub = solve(op.sc_op) 
ui = MultiField.backward_static_condensation(op,ub)
