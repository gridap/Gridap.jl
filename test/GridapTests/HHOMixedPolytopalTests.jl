module HHOMixedPolytopalTests

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

function reconstruction_operator(ptopo,order,X,Ω,Γp,dΩp,dΓp)
  L = FESpaces.PolytopalFESpace(Ω, Float64, order+1; space=:P)
  Λ = FESpaces.PolytopalFESpace(Ω, Float64, 0; space=:P)

  n = get_normal_vector(Γp)
  Πn(v) = ∇(v)⋅n
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((u,λ),(v,μ))   = ∫( (∇(u)⋅∇(v)) + (μ*u) + (λ*v) )dΩp
  rhs((uT,uF),(v,μ)) =  ∫( (∇(uT)⋅∇(v)) + (uT*μ) )dΩp + ∫( (uF - Π(uT,Γp))*(Πn(v)) )dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  mfs = Y.multi_field_style
  W = MultiFieldFESpace([L,Λ];style=mfs)
  R = LocalOperator(
    LocalPenaltySolveMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return R
end


##############################################################
#u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
u(x) = x[1] + x[2]
f(x) = -Δ(u)(x)

D = 3
nc = Tuple(fill(2,D))
domain = Tuple(repeat([0,1],D))
model = simplexify(CartesianDiscreteModel(domain,nc);positive=true)

if D == 2
  vmodel = Gridap.Geometry.voronoi(model)
else
  vmodel = Geometry.PolytopalDiscreteModel(model)
end

Ω = Triangulation(ReferenceFE{D}, vmodel)
Γ = Triangulation(ReferenceFE{D-1}, vmodel)

ptopo = Geometry.PatchTopology(vmodel)
Ωp = Geometry.PatchTriangulation(vmodel,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(vmodel,ptopo)

order = 1
qdegree = 2*(order+1)

dΩp = Measure(Ωp,qdegree)
dΓp = Measure(Γp,qdegree)

V = FESpaces.PolytopalFESpace(Ω, Float64, order+1; space=:P) # Bulk space
if D == 2
  reffe_M = ReferenceFE(lagrangian, Float64, order; space=:P) 
  M = FESpace(Γ, reffe_M; conformity=:L2, dirichlet_tags="boundary") # Skeleton space
else
  M = FESpaces.PolytopalFESpace(Γ, Float64, order; space=:P, dirichlet_tags="boundary")
end
N = TrialFESpace(M,u)

mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
X   = MultiFieldFESpace([V, N];style=mfs)
Y   = MultiFieldFESpace([V, M];style=mfs) 
Xp  = FESpaces.PatchFESpace(X,ptopo)

PΓ = projection_operator(M, Γp, dΓp)
R  = reconstruction_operator(ptopo,order,Y,Ωp,Γp,dΩp,dΓp)

global_assem = SparseMatrixAssembler(X,Y)
patch_assem = FESpaces.PatchAssembler(ptopo,X,Y)

function a(u,v)
  Ru_Ω, Ru_Γ = R(u)
  Rv_Ω, Rv_Γ = R(v)
  return ∫(∇(Ru_Ω)⋅∇(Rv_Ω) + ∇(Ru_Γ)⋅∇(Rv_Ω) + ∇(Ru_Ω)⋅∇(Rv_Γ) + ∇(Ru_Γ)⋅∇(Rv_Γ))dΩp
end

hT = collect(get_array(∫(1)dΩp))
hTinv = CellField(1 ./ hT,Ωp)
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

# Monolithic solve
A, b = weakform()
x = A \ b

ui, ub = FEFunction(X,x)
eu  = ui - u 
l2u = sqrt(sum( ∫(eu * eu)dΩp))
h1u = l2u + sqrt(sum( ∫(∇(eu) ⋅ ∇(eu))dΩp))
@test l2u < 1e-10

# Static condensation
op = MultiField.StaticCondensationOperator(X,patch_assem,patch_weakform())
ui, ub = solve(op) 

eu  = ui - u
l2u = sqrt(sum( ∫(eu * eu)dΩp))
h1u = l2u + sqrt(sum( ∫(∇(eu) ⋅ ∇(eu))dΩp))
@test l2u < 1e-10

end