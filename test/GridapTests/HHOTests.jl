module HHOTests

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

function mf_projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((p,),(q,)) = ∫(p⋅q)dΩ
  rhs((uT,uF),(q,)) = ∫(Π(uT,Ω)⋅q + Π(uF,Ω)⋅q)dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  W = MultiFieldFESpace([V0],style=MultiField.BlockMultiFieldStyle())
  P = LocalOperator(
    LocalSolveMap(), W, lhs, rhs; trian_out = Ω, space_out = V
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
  P = LocalOperator(
    LocalSolveMap(), ptopo, W, Y, lhs, rhs; trian_out = Ω, space_out = V
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
  R = LocalOperator(
    LocalPenaltySolveMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return R
end

##############################################################
# uex(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
uex(x) = x[1]+x[2]
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
PΩ_mf = patch_projection_operator(ptopo,V,Yp,Ωp,dΩp)

global_assem = SparseMatrixAssembler(X,Y)
patch_assem = FESpaces.PatchAssembler(ptopo,X,Y)

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
  PΓRu_Ω, PΓRu_Γ = PΓ_mf(Ru)
  PΓPΩRu_Ω, PΓPΩRu_Γ = PΓ_mf(PΩ_mf(Ru))
  return (PΓRu_Ω - PΓPΩRu_Ω) + (PΓRu_Γ - PΓPΩRu_Γ)
end

function weakform()
  u, v = get_trial_fe_basis(X), get_fe_basis(X)

  mass_Γ(u,v) = ∫(hFinv*(u*v) )dΓp
  Ru, Rv = R(u), R(v)
  SΓa_u, SΓa_v = SΓa(u), SΓa(v)
  SΓb_u, SΓb_v = SΓb(Ru), SΓb(Rv)

  data = FESpaces.collect_and_merge_cell_matrix_and_vector(
    (Xp, Xp, a(Ru,Rv) + mass_Γ(SΓb_u,SΓb_v), DomainContribution(), zero(Xp)),
    (X, X, mass_Γ(SΓa_u,SΓa_v), l(v), zero(X)),
    (X, Xp, mass_Γ(SΓa_u,SΓb_v), DomainContribution(), zero(X)),
    (Xp, X, mass_Γ(SΓb_u,SΓa_v), DomainContribution(), zero(Xp)),
  )

  return assemble_matrix_and_vector(global_assem,data)
end

function patch_weakform()
  u, v = get_trial_fe_basis(X), get_fe_basis(X);

  mass_Γ(u,v) = ∫(hFinv*(u*v))dΓp
  Ru, Rv = R(u), R(v)
  SΓa_u, SΓa_v = SΓa(u), SΓa(v)
  SΓb_u, SΓb_v = SΓb(Ru), SΓb(Rv)

  data = FESpaces.collect_and_merge_cell_matrix_and_vector(patch_assem,
    (Xp, Xp, a(Ru,Rv) + mass_Γ(SΓb_u,SΓb_v), DomainContribution(), zero(Xp)),
    (X, X, mass_Γ(SΓa_u,SΓa_v), l(v), zero(X)),
    (X, Xp, mass_Γ(SΓa_u,SΓb_v), DomainContribution(), zero(X)),
    (Xp, X, mass_Γ(SΓb_u,SΓa_v), DomainContribution(), zero(Xp)),
  )

  return assemble_matrix_and_vector(patch_assem,data)
end

A, b = weakform()
x = A \ b
ui, ub = FEFunction(X,x);

eh = ui - uex
el2 = sqrt(sum(∫(eh⋅eh)dΩp))
@test el2 < 1e-10

# Static condensation
op = MultiField.StaticCondensationOperator(X,patch_assem,patch_weakform())
xh = solve(op)
uΩ, uΓ = xh

eu = uΩ - uex 
el2 = sqrt(sum( ∫(eu * eu)dΩp))
@test el2 < 1e-10

Rxh = R(xh)
eu = Rxh - uex
el2 = sqrt(sum( ∫(eu * eu)dΩp))
@test el2 < 1e-10

end