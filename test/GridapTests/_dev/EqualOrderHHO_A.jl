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

##############################################################
uex(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
f(x) = -Δ(u)(x)

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
PΩ = projection_operator(V, Ωp, dΩp)

u = get_trial_fe_basis(X);
v = get_fe_basis(Y);

# SΓa
function SΓa(u)
  u_Ω, u_Γ = u
  return PΓ(u_Ω) - u_Γ
end

uT, uF = u;
@enter PΓ(uT)
Sa = SΓa(u)

@enter Ru = R(u)

Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
mass(u,v,Ω,dΩ) = ∫(u⋅Π(v,Ω))dΩ
V0 = FESpaces.FESpaceWithoutBCs(V)

pu = get_trial_fe_basis(V0)
pv = get_fe_basis(V)
lhs = mass(pu, pv, Ωp, dΩp)
RuΩ, RuΓ = Ru
# rhs = mass(RuΩ, pv, Ωp, dΩp) + mass(RuΓ, pv, Ωp, dΩp)
rhs = ∫( (RuΩ+RuΓ) * pv)dΩp 

_RuΩ = FESpaces.SingleFieldFEBasis(lazy_map(transpose,RuΩ.cell_field.args[1].args[1]),Ωp,FESpaces.TrialBasis(),ReferenceDomain())
_RuΓ = FESpaces.SingleFieldFEBasis(lazy_map(transpose,RuΓ.cell_field.args[1].args[1]),Ωp,FESpaces.TrialBasis(),ReferenceDomain())

rhs_Ω = ∫(_RuΩ* pv)dΩp
rhs_Γ = ∫(_RuΓ* pv)dΩp
get_array(rhs_Γ)[1]

mean_Γ = get_array(∫(_RuΓ)dΩp)

lhs_assem = FESpaces.PatchAssembler(ptopo,V0,V0)
lhs_mats  = FESpaces.assemble_matrix(
  lhs_assem,FESpaces.collect_patch_cell_matrix(lhs_assem,V0,V0,lhs)
)

rhs_assem_Ω = FESpaces.PatchAssembler(ptopo,V,V0)
rhs_mats_Ω  = FESpaces.assemble_matrix(
  rhs_assem_Ω,FESpaces.collect_patch_cell_matrix(rhs_assem_Ω,V,V0,rhs_Ω)
)

M0 = FESpaces.FESpaceWithoutBCs(M)
M0p = FESpaces.PatchFESpace(M0,ptopo)
rhs_assem_Γ = FESpaces.PatchAssembler(ptopo,M0p,V0)
rhs_mats_Γ = FESpaces.assemble_matrix(
  rhs_assem_Γ,FESpaces.collect_patch_cell_matrix(rhs_assem_Γ,M0p,V0,rhs_Γ)
)

mass_Ω(u,v) = mass(u,v,Ωp,dΩp)
PΩRΩ = FESpaces.LocalOperator(FESpaces.LocalSolveMap(), ptopo, V0, V, mass_Ω, mass_Ω)
PΩRΓ = FESpaces.LocalOperator(FESpaces.LocalSolveMap(), ptopo, V0, M0p, mass_Ω, mass_Ω)

PΩRuΩ = PΩRΩ(_RuΩ)
PΩRuΓ = PΩRΓ(_RuΓ)

wΩ = FESpaces.SingleFieldFEBasis(CellData.get_data(_RuΩ - PRuΩ),Ωp,FESpaces.TrialBasis(),ReferenceDomain())
wΓ = FESpaces.SingleFieldFEBasis(CellData.get_data(_RuΓ - PRuΓ),Ωp,FESpaces.TrialBasis(),ReferenceDomain())

PΓwΩ = PΓ(wΩ)
PΓwΓ = PΓ(wΓ)

_PΓwΩ = FESpaces.SingleFieldFEBasis(CellData.get_data(PΓwΩ),Γp,FESpaces.TrialBasis(),ReferenceDomain())
_PΓwΓ = FESpaces.SingleFieldFEBasis(CellData.get_data(PΓwΓ),Γp,FESpaces.TrialBasis(),ReferenceDomain())

mf_PΓwΩ = MultiField.MultiFieldFEBasisComponent(_PΓwΩ,1,2)
mf_PΓwΓ = MultiField.MultiFieldFEBasisComponent(_PΓwΓ,2,2)


uΩ, uΓ = u;
PΓuΩ = PΓ(uΩ)  
∫(PΓuΩ + mf_PΓwΩ + mf_PΓwΓ)dΓp


# P = FESpaces.LocalOperator(
#   FESpaces.LocalSolveMap(), V0, mass, mass; trian_out = Ω
# )











global_assem = SparseMatrixAssembler(X,Y)
patch_assem = FESpaces.PatchAssembler(ptopo,X,Y)

function a(u,v)
  Ru_Ω, Ru_Γ = R(u)
  Rv_Ω, Rv_Γ = R(v)
  return ∫(∇(Ru_Ω)⋅∇(Rv_Ω) + ∇(Ru_Γ)⋅∇(Rv_Ω) + ∇(Ru_Ω)⋅∇(Rv_Γ) + ∇(Ru_Γ)⋅∇(Rv_Γ))dΩp
end

hTinv =  CellField(1 ./ (sqrt(2).*sqrt.(get_array(∫(1)dΩp))),Ωp)
function s(u,v)
  function SΓ(u)
    u_Ω, u_Γ = u
    return PΓ(u_Ω) - u_Γ
  end
  return ∫(hTinv * (SΓ(u)⋅SΓ(v)))dΓp
end

l((vΩ,vΓ)) = ∫(f⋅vΩ)dΩp

S = ∫(hTinv*SΓa(u))*dΓp

data2 = FESpaces.collect_cell_matrix_and_vector(Xp,Xp,a(u,v),DomainContribution(),zero(Xp))

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
