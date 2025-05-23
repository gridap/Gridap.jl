
"""
References: 

- A hybrid high-order locking-free method for linear elasticity on general meshes. Di Pietro et al. 2014
- Locking-free hybrid high-order method for linear elasticity. Carstensen et al. 2024

TODO: 

- I believe the papers use regular HHO. Can we find a reference for mixed order? 
- Check convergence rates
"""

using Test
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs, Gridap.TensorValues

skew_symmetric_gradient(f) = Operation(TensorValues.skew_symmetric_part)(gradient(f))
ν(f) = skew_symmetric_gradient(f)

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
  reffe_Λ1 = ReferenceFE(lagrangian, VectorValue{2,Float64}, 0)
  reffe_Λ2 = ReferenceFE(lagrangian, Float64, 0) # In 3D, should be type VectorValue{3,Float64}
  Λ1 = FESpace(Ω, reffe_Λ1; conformity=:L2)
  Λ2 = FESpace(Ω, reffe_Λ2; conformity=:L2)

  nrel = get_normal_vector(Γp)
  σ(ν_u,λ) = ν_u[1,2]*λ # Constraint on the off-diagonal term
  Πn(v) = nrel⋅∇(v)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((u,λ1,λ2),(v,μ1,μ2)) = ∫((ε(u)⊙ε(v)) + (μ1⋅u) + (λ1⋅v) + σ∘(ν(u),μ2) + σ∘(ν(v),λ2))dΩp
  rhs((uT,uF),(v,μ1,μ2)) = ∫((ε(uT)⊙ε(v)) + (uT⋅μ1))dΩp + ∫((uF - Π(uT,Γp))⋅(Πn(v)) + σ∘(outer(uF,Π(uT,Γp)),μ2))dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  mfs = MultiField.BlockMultiFieldStyle(2,(1,2))
  W = MultiFieldFESpace([L,Λ1,Λ2];style=mfs)
  R = FESpaces.LocalOperator(
    FESpaces.HHO_ReconstructionOperatorMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return R
end

##############################################################
# u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
u(x) = VectorValue(sin(2*π*x[1])*sin(2*π*x[2]),cos(2*π*x[1])*cos(2*π*x[2]))
f(x) = -(∇⋅(ε(u)))(x)

D = 2
n = 64
model = UnstructuredDiscreteModel(CartesianDiscreteModel(Tuple(repeat([0,1],D)),ntuple(i -> n, D)))
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

global_assem = SparseMatrixAssembler(X,Y)
patch_assem = FESpaces.PatchAssembler(ptopo,X,Y)

function a(u,v)
  Ru_Ω, Ru_Γ = R(u)
  Rv_Ω, Rv_Γ = R(v)
  return ∫(ε(Ru_Ω)⊙ε(Rv_Ω) + ε(Ru_Γ)⊙ε(Rv_Ω) + ε(Ru_Ω)⊙ε(Rv_Γ) + ε(Ru_Γ)⊙ε(Rv_Γ))dΩp
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
# A, b = weakform()
# x = A \ b

# Static condensation
op = MultiField.StaticCondensationOperator(X,V,N,patch_assem,patch_weakform())

ub = solve(op.sc_op) 
ui = MultiField.backward_static_condensation(op,ub)

dΩ = Measure(Ω,qdegree)
l2_ui = sqrt(sum(∫((ui - u)⋅(ui - u))*dΩ))
