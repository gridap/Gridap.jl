
"""
References: 

- A hybrid high-order locking-free method for linear elasticity on general meshes. Di Pietro et al. 2014
- Locking-free hybrid high-order method for linear elasticity. Carstensen et al. 2024
- The hybrid high-order method for polytopal meshes. Di Pietro, Droniou. 2020

"""

using Test
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs, Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields: skew_symmetric_gradient

ν(f) = skew_symmetric_gradient(f)

function projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  P = LocalOperator(
    LocalSolveMap(), V0, mass, mass; trian_out = Ω
  )
  return P
end

function reconstruction_operator(ptopo,L,X,Ω,Γp,dΩp,dΓp)
  reffe_Λ1 = ReferenceFE(lagrangian, VectorValue{2,Float64}, 0)
  reffe_Λ2 = ReferenceFE(lagrangian, VectorValue{3,Float64}, 0)
  Λ1 = FESpace(Ω, reffe_Λ1; conformity=:L2)
  Λ2 = FESpace(Ω, reffe_Λ2; conformity=:L2)

  nrel = get_normal_vector(Γp)
  
  Πn(v) = nrel⋅ε(v)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))

  σ(ν_u,λ) = ν_u[1,1]*λ[1] + ν_u[1,2]*λ[2] + ν_u[2,2]*λ[3] 
  β(a,b) = -0.5*(a⊗b - b⊗a)

  lhs((u,λ1,λ2),(v,μ1,μ2)) = ∫((ε(u)⊙ε(v)) + (μ1⋅u) + (λ1⋅v) + σ∘(ν(u),μ2) + σ∘(ν(v),λ2))dΩp
  rhs((uT,uF),(v,μ1,μ2)) = ∫((ε(uT)⊙ε(v)) + (uT⋅μ1))dΩp + ∫((uF - Π(uT,Γp))⋅(Πn(v)) + σ∘(β∘(uF,nrel),μ2))dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  mfs = MultiField.BlockMultiFieldStyle(2,(1,2))
  W = MultiFieldFESpace([L,Λ1,Λ2];style=mfs)
  R = LocalOperator(
    LocalPenaltySolveMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return R
end

function divergence_operator(ptopo,L,X,Ω,Γp,dΩp,dΓp)
  nrel = get_normal_vector(Γp)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((p,),(q,))   =  ∫(p⋅q)dΩp
  rhs((uT,uF),(q,)) =  ∫((∇⋅uT)⋅q)dΩp + ∫((uF - Π(uT,Γp))⋅nrel * q)dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  mfs = MultiField.BlockMultiFieldStyle(1)
  W = MultiFieldFESpace([L];style=mfs)
  D = LocalOperator(
    LocalSolveMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return D
end

function gradient_operator(ptopo,X,order,Ω,Γp,dΩp,dΓp)
  reffe_L = ReferenceFE(lagrangian,SymTensorValue{2,Float64}, order; space = :P)
  L = FESpace(Ω, reffe_L; conformity=:L2)

  nΓ = get_normal_vector(Γp)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((p,),(q,))   =  ∫(p⊙q)dΩp
  rhs((uT,uF),(q,)) =  ∫(ε(uT)⊙q)dΩp + ∫((uF - Π(uT,Γp)) ⋅ (q⋅nΓ))dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  W = MultiFieldFESpace([L];style=MultiField.BlockMultiFieldStyle(1))
  G = LocalOperator(
    LocalSolveMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return G
end

##############################################################
# u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])

μ = 1.0
λ = 1.0
I2 = one(SymTensorValue{2,Float64})
I4 = one(SymFourthOrderTensorValue{2,Float64})
C = λ * I2⊗I2 + μ * I4 # Isotropic elasticity tensor

u(x) = VectorValue(sin(2*π*x[1])*sin(2*π*x[2]),cos(2*π*x[1])*cos(2*π*x[2]))
g(x) = C⊙(ε(u)(x))
f(x) = -(∇⋅g)(x)


D = 2
n = 64
model = UnstructuredDiscreteModel(CartesianDiscreteModel(Tuple(repeat([0,1],D)),ntuple(i -> n, D)))
D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D}, model)
Γ = Triangulation(ReferenceFE{D-1}, model)

ptopo = Geometry.PatchTopology(model)
Ωp = Geometry.PatchTriangulation(model,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)

order = 2
qdegree = 2*(order+1)

dΩp = Measure(Ωp,qdegree)
dΓp = Measure(Γp,qdegree)

##########################
# Mixed order variant
##########################
reffe_V = ReferenceFE(lagrangian, VectorValue{D,Float64}, order+1; space=:P)   # Bulk space
reffe_M = ReferenceFE(lagrangian, VectorValue{D,Float64}, order; space=:P)     # Skeleton space
reffe_L = ReferenceFE(lagrangian, VectorValue{D,Float64}, order+1; space=:P)   # Reconstruction space
V = FESpace(Ω, reffe_V; conformity=:L2)
M = FESpace(Γ, reffe_M; conformity=:L2, dirichlet_tags="boundary")
L = FESpace(Ω, reffe_L; conformity=:L2)
N = TrialFESpace(M,u)

mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
X  = MultiFieldFESpace([V, N];style=mfs)
Y  = MultiFieldFESpace([V, M];style=mfs)
Xp = FESpaces.PatchFESpace(X,ptopo)

PΓ = projection_operator(M, Γp, dΓp)
R  = reconstruction_operator(ptopo,L,Y,Ωp,Γp,dΩp,dΓp)
G  = gradient_operator(ptopo,X,order,Ωp,Γp,dΩp,dΓp)

global_assem = SparseMatrixAssembler(X,Y)
patch_assem = FESpaces.PatchAssembler(ptopo,X,Y)

function a(u,v)
  Gu_Ω, Gu_Γ = G(u)
  Gv_Ω, Gv_Γ = G(v)
  Gu = Gu_Ω + Gu_Γ
  Gv = Gv_Ω + Gv_Γ
  return ∫(Gu⊙C⊙Gv)dΩp
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

function weakform()
  u, v = get_trial_fe_basis(X), get_fe_basis(Y)
  data = FESpaces.collect_and_merge_cell_matrix_and_vector(
    (Xp, Xp, a(u,v), DomainContribution(), zero(Xp)),
    (X, Y, s(u,v), l(v), zero(X))
  )
  assemble_matrix_and_vector(global_assem,data)
end

function patch_weakform()
  u, v = get_trial_fe_basis(X), get_fe_basis(Y)
  data = FESpaces.collect_and_merge_cell_matrix_and_vector(patch_assem,
    (Xp, Xp, a(u,v), DomainContribution(), zero(Xp)),
    (X, Y, s(u,v), l(v), zero(X))
  )
  return assemble_matrix_and_vector(patch_assem,data)
end

# Static condensation
op = MultiField.StaticCondensationOperator(X,patch_assem,patch_weakform())
xh = solve(op)
ui, ub = xh

dΩ = Measure(Ω,qdegree)
l2_ui = sqrt(sum(∫((ui - u)⋅(ui - u))*dΩ))

Ru = R(xh)
dΩ = Measure(Ω,qdegree)
l2_ui = sqrt(sum(∫((Ru - u)⋅(Ru - u))*dΩ))
