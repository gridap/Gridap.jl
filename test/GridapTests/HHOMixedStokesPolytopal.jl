"""
References: 

- A discontinuous skeletal method for the viscosity-dependent Stokes problem. Di Pietro et al. 2016
- The hybrid high-order method for polytopal meshes. Di Pietro, Droniou. 2020

"""

using Test
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs
using Gridap.Arrays

function l2_error(uh,u,dΩ)
  eh = uh - u
  return sqrt(sum(∫(eh⋅eh)*dΩ))
end

function swap_field_ids(u, ids, nfields)
  @assert length(ids) == length(u)
  fields = map((ui,id) -> swap_field_ids(ui,id,nfields), u, ids)
  return MultiField.MultiFieldCellField(fields)
end

function swap_field_ids(u::MultiField.MultiFieldFEBasisComponent, id, nfields)
  return MultiField.MultiFieldFEBasisComponent(u.single_field, id, nfields)
end

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
  D = num_cell_dims(Ω)
  Λ = FESpaces.PolytopalFESpace(Ω, VectorValue{D,Float64}, 0)

  nrel = get_normal_vector(Γp)
  Πn(v) = nrel⋅∇(v)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((u,λ),(v,μ))   =  ∫( (∇(u)⊙∇(v)) + (μ⋅u) + (λ⋅v) )dΩp
  rhs((uT,uF),(v,μ)) =  ∫( (∇(uT)⊙∇(v)) + (uT⋅μ) )dΩp + ∫( (uF - Π(uT,Γp))⋅(Πn(v)) )dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
  W = MultiFieldFESpace([L,Λ];style=mfs)
  R = LocalOperator(
    LocalPenaltySolveMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  _R(u) = swap_field_ids(R(swap_field_ids(u,[1,2],2)),[1,3],5)
  return _R
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
  _D(u) = swap_field_ids(D(swap_field_ids(u,[1,2],2)),[1,3],5)
  return _D
end

##############################################################
# u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
u(x) = VectorValue(x[1],-x[2])
p(x) = x[1] + x[2]
f(x) = -Δ(u)(x) - ∇(p)(x)

n = 8
model = Geometry.voronoi(simplexify(CartesianDiscreteModel((0,1,0,1),(n,n))))

D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D}, model)
Γ = Triangulation(ReferenceFE{D-1}, model)

ptopo = Geometry.PatchTopology(model)
Ωp = Geometry.PatchTriangulation(model,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)

order = 1
qdegree = 2*(order+1)

dΩ = Measure(Ω,qdegree)
dΩp = Measure(Ωp,qdegree)
dΓp = Measure(Γp,qdegree)

##########################
# Mixed order variant
##########################

V = FESpaces.PolytopalFESpace(Ω, VectorValue{2,Float64}, order+1)
M = FESpaces.PolytopalFESpace(Γ, VectorValue{2,Float64}, order; dirichlet_tags="boundary")
L = FESpaces.PolytopalFESpace(Ω, VectorValue{2,Float64}, order+1)
Q = FESpaces.PolytopalFESpace(Ω, Float64, order) # Pressure space

Q̂ = FESpaces.PolytopalFESpace(Ω, Float64, order; local_kernel=:constants)
Q̄ = FESpaces.PolytopalFESpace(Ω, Float64, 0) 

Λ = ConstantFESpace(Ω) # Lagrange multiplier for zero mean
N = TrialFESpace(M,u)

mfs = MultiField.BlockMultiFieldStyle(2,(2,3))
X  = MultiFieldFESpace([V, Q̂, N, Q̄, Λ];style=mfs)
Y  = MultiFieldFESpace([V, Q̂, M, Q̄, Λ];style=mfs)
Xp = FESpaces.PatchFESpace(X,ptopo)

YR = MultiFieldFESpace([V, M];style=MultiField.BlockMultiFieldStyle(2,(1,1)))

PΓ = projection_operator(M, Γp, dΓp)
R  = reconstruction_operator(ptopo,L,YR,Ωp,Γp,dΩp,dΓp)
D  = divergence_operator(ptopo,Q,YR,Ωp,Γp,dΩp,dΓp)

global_assem = SparseMatrixAssembler(X,Y)
patch_assem = FESpaces.PatchAssembler(ptopo,X,Y)

function a(Ru,Rv)
  u_Ω, u_Γ = Ru
  v_Ω, v_Γ = Rv
  return ∫(∇(u_Ω)⊙∇(v_Ω) + ∇(u_Γ)⊙∇(v_Ω) + ∇(u_Ω)⊙∇(v_Γ) + ∇(u_Γ)⊙∇(v_Γ))dΩp
end

function b(Du,q)
  Du_Ω, Du_Γ = Du
  return ∫((Du_Ω + Du_Γ)⋅q)dΩp
end

function c(q,λ)
  return ∫(q⋅λ)dΩp
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

function weakform(x,y)
  u, p̂, p̄, λ = (x[1],x[3]), x[2], x[4], x[5]
  v, q̂, q̄, μ = (y[1],y[3]), y[2], y[4], y[5]
  p, q = p̂ + p̄, q̂ + q̄

  Ru, Rv = R(u), R(v)
  Du, Dv = D(u), D(v)

  data = FESpaces.collect_and_merge_cell_matrix_and_vector(
    (Xp, Xp, a(Ru,Rv), DomainContribution(), zero(Xp)),
    (X, Y, s(u,v)+c(p̄,μ)+c(q̄,λ), l(v), zero(X)),
    (X, Xp, b(Dv,p), DomainContribution(), zero(X)),
    (Xp, X, b(Du,q), DomainContribution(), zero(Xp))
  )
  assemble_matrix_and_vector(global_assem,data)
end

function patch_weakform(x,y)
  u, p̂, p̄, λ = (x[1],x[3]), x[2], x[4], x[5]
  v, q̂, q̄, μ = (y[1],y[3]), y[2], y[4], y[5]
  p, q = p̂ + p̄, q̂ + q̄

  Ru, Rv = R(u), R(v)
  Du, Dv = D(u), D(v)

  data = FESpaces.collect_and_merge_cell_matrix_and_vector(patch_assem,
    (Xp, Xp, a(Ru,Rv), DomainContribution(), zero(Xp)),
    (X, Y, s(u,v)+c(p̄,μ)+c(q̄,λ), l(v), zero(X)),
    (X, Xp, b(Dv,p), DomainContribution(), zero(X)),
    (Xp, X, b(Du,q), DomainContribution(), zero(Xp))
  )
  assemble_matrix_and_vector(patch_assem,data)
end

x, y = get_trial_fe_basis(X), get_fe_basis(Y)

# Monolithic solve
A, B = weakform(x,y)
xh = FEFunction(X, A \ B)
uhT, uhF, ph = xh[1], xh[3], xh[2]+xh[4]

l2_uT = l2_error(uhT,u,dΩp)
l2_uF = l2_error(uhF,u,dΓp)
l2_p = l2_error(ph,p,dΩp) # Because p is not zero mean

# Static condensation
op = MultiField.StaticCondensationOperator(
  X,patch_assem,patch_weakform(x,y)
)

xh = solve(op)
ui, ub, ph = xh[1], xh[3], xh[2]+xh[4]

l2_ui = sqrt(sum(∫((ui - u)⋅(ui - u))*dΩ))
@test l2_ui < 1e-10
