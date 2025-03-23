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

function mf_projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((p,),(q,)) = ∫(p⋅q)dΩ
  rhs((uT,uF),(q,)) = ∫(Π(uT,Ω)⋅q + Π(uF,Ω)⋅q)dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  W = MultiFieldFESpace([V0],style=MultiField.BlockMultiFieldStyle())
  P = FESpaces.LocalOperator(
    FESpaces.LocalSolveMap(), W, lhs, rhs; trian_out = Ω, space_out = V
  )
  return P
end

function patch_projection_operator(ptopo,V,X,Ω,dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((p,),(q,)) = ∫(p⋅q)dΩ
  rhs((uT,uF),(q,)) = ∫(Π(uT,Ω)⋅q + Π(uF,Ω)⋅q)dΩ

  V0 = FESpaces.FESpaceWithoutBCs(V)
  Y = FESpaces.FESpaceWithoutBCs(X)
  W = MultiFieldFESpace([V0],style=MultiField.BlockMultiFieldStyle())
  P = FESpaces.LocalOperator(
    FESpaces.LocalSolveMap(), ptopo, W, Y, lhs, rhs; trian_out = Ω, space_out = V
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
  W = MultiFieldFESpace([L,Λ];style=MultiField.MultiFieldStyle(Y))
  R = FESpaces.LocalOperator(
    FESpaces.HHO_ReconstructionOperatorMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return R
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

order = 0
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
Yp = FESpaces.PatchFESpace(Y,ptopo)

R  = reconstruction_operator(ptopo,L,Y,Ωp,Γp,dΩp,dΓp)
PΓ = projection_operator(M, Γp, dΓp)
PΓ_mf = mf_projection_operator(M, Γp, dΓp)
PΩ_mf = patch_projection_operator(ptopo,V,Xp,Ωp,dΩp)

global_assem = SparseMatrixAssembler(X,Y)
patch_assem = FESpaces.PatchAssembler(ptopo,X,Y)

u, v = get_trial_fe_basis(X), get_fe_basis(Y);

Ru, Rv = R(u), R(v)

function a(Ru,Rv)
  Ru_Ω, Ru_Γ = Ru
  Rv_Ω, Rv_Γ = Rv
  return ∫(∇(Ru_Ω)⋅∇(Rv_Ω) + ∇(Ru_Γ)⋅∇(Rv_Ω) + ∇(Ru_Ω)⋅∇(Rv_Γ) + ∇(Ru_Γ)⋅∇(Rv_Γ))dΩp
end

l((vΩ,vΓ)) = ∫(f⋅vΩ)dΩp

hFinv =  CellField(1 ./ get_array(∫(1)dΓp) ,Γp)
function SΓa(u)
  u_Ω, u_Γ = u
  return PΓ(u_Ω) - u_Γ 
end
function SΓb(Ru)
  function minus(a::MultiField.MultiFieldFEBasisComponent,b::MultiField.MultiFieldFEBasisComponent)
    sf = a.single_field - b.single_field
    sf_basis = FESpaces.SingleFieldFEBasis(CellData.get_data(sf),get_triangulation(sf),BasisStyle(a),DomainStyle(sf))
    MultiField.MultiFieldFEBasisComponent(sf_basis,a.fieldid,a.nfields)
  end
  minus(a,b) = MultiFieldCellField(map(minus,a,b))
  SΓb_Ω, SΓb_Γ = PΓ_mf(minus(Ru,PΩ_mf(Ru)))
  return SΓb_Ω + SΓb_Γ
end

function weakform(u,v)
  mass_Γ(u,v) = ∫(hFinv*(u*v) )dΓp
  Ru, Rv = R(u), R(v)
  SΓa_u, SΓa_v = SΓa(u), SΓa(v)
  SΓb_u, SΓb_v = SΓb(Ru), SΓb(Rv)
  ll = FESpaces.collect_cell_matrix_and_vector(Xp,Xp,a(Ru,Rv),DomainContribution(),zero(Xp))
  aa = FESpaces.collect_cell_matrix_and_vector(X,X,mass_Γ(SΓa_u,SΓa_v),l(v),zero(X))
  ab = FESpaces.collect_cell_matrix_and_vector(X,Xp,mass_Γ(SΓa_u,SΓb_v),DomainContribution(),zero(X))
  ba = FESpaces.collect_cell_matrix_and_vector(Xp,X,mass_Γ(SΓb_u,SΓa_v),DomainContribution(),zero(Xp))
  bb = FESpaces.collect_cell_matrix_and_vector(Xp,Xp,mass_Γ(SΓb_u,SΓb_v),DomainContribution(),zero(Xp))
  data = FESpaces.merge_assembly_matvec_data(ll,aa,ab,ba,bb)
  return assemble_matrix_and_vector(global_assem,data)
end

A, b = weakform(u,v)
x = A \ b

ub, us = FEFunction(X,x);

eh = ub - uex
sqrt(sum(∫(eh⋅eh)dΩp))

function patch_weakform(u,v)
  mass_Γ(u,v) = ∫(hFinv*(u*v))dΓp
  Ru, Rv = R(u), R(v)
  SΓa_u, SΓa_v = SΓa(u), SΓa(v)
  SΓb_u, SΓb_v = SΓb(Ru), SΓb(Rv)
  ll = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,Xp,Xp,a(Ru,Rv),DomainContribution(),zero(Xp))
  aa = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,X,X,mass_Γ(SΓa_u,SΓa_v),l(v),zero(X))
  ab = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,X,Xp,mass_Γ(SΓa_u,SΓb_v),DomainContribution(),zero(X))
  ba = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,Xp,X,mass_Γ(SΓb_u,SΓa_v),DomainContribution(),zero(Xp))
  bb = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,Xp,Xp,mass_Γ(SΓb_u,SΓb_v),DomainContribution(),zero(Xp))
  data = FESpaces.merge_assembly_matvec_data(ll,aa,ab,ba,bb)
  return assemble_matrix_and_vector(patch_assem,data)
end

# Static condensation
op = MultiField.StaticCondensationOperator(X,V,N,patch_assem,patch_weakform(u,v))

uΓ = solve(op.sc_op) 
uΩ = MultiField.backward_static_condensation(op,uΓ)

eu = uΩ - uex 
sqrt(sum( ∫(eu * eu)dΩp))
