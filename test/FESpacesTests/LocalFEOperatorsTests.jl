module LocalFEOperatorsTests

# A `LocalOperator` applied to a space `X` returns the space `R(X)`: the dofs of `X`, with
# the images of the basis of `X` as shape functions.

using Test
using LinearAlgebra
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Arrays

Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))

# Cell-wise, single-field: L2 projection onto V
function projection_operator(V, Ω, dΩ)
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  return LocalOperator(LocalSolveMap(), V0, mass, mass; trian_out = Ω)
end

# Cell-wise, multi-field input
function mf_projection_operator(V, Ω, dΩ)
  lhs((p,),(q,)) = ∫(p⋅q)dΩ
  rhs((uT,uF),(q,)) = ∫(Π(uT,Ω)⋅q + Π(uF,Ω)⋅q)dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  W = MultiFieldFESpace([V0],style=MultiField.BlockMultiFieldStyle())
  return LocalOperator(LocalStagedSolveMap((((1,),(1,)),)), W, lhs, rhs; trian_out = Ω, space_out = V)
end

# Patch-wise, multi-field input: HHO reconstruction
function reconstruction_operator(ptopo,L,X,Ω,Γp,dΩp,dΓp)
  Λ = FESpaces.PolytopalFESpace(Ω, Float64, 0; space=:P)
  nrel = get_normal_vector(Γp)
  Πn(v) = ∇(v)⋅nrel
  lhs((u,λ),(v,μ))   =  ∫( (∇(u)⋅∇(v)) + (μ*u) + (λ*v) )dΩp
  rhs((uT,uF),(v,μ)) =  ∫( (∇(uT)⋅∇(v)) + (uT*μ) )dΩp + ∫( (uF - Π(uT,Γp))*(Πn(v)) )dΓp
  Y = FESpaces.FESpaceWithoutBCs(X)
  W = MultiFieldFESpace([L,Λ];style=Y.multi_field_style)
  return LocalOperator(LocalPenaltySolveMap(), ptopo, W, Y, lhs, rhs; space_out = L)
end

l2(e,dΩ) = sqrt(sum(∫(e⋅e)dΩ))

g(x) = 1 + x[1] + sin(3*x[2]) # Non-homogeneous Dirichlet data
order = 1

model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(3,3)))
vmodel = Geometry.voronoi(simplexify(model))
D = num_cell_dims(vmodel)
Ω = Triangulation(ReferenceFE{D}, vmodel)
Γ = Triangulation(ReferenceFE{D-1}, vmodel)

ptopo = Geometry.PatchTopology(vmodel)
Ωp = Geometry.PatchTriangulation(vmodel,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(vmodel,ptopo)

qdegree = 2*(order+1)
dΩp = Measure(Ωp,qdegree)
dΓp = Measure(Γp,qdegree)

reffe_M = ReferenceFE(lagrangian, Float64, order; space=:P)
V = FESpaces.PolytopalFESpace(Ω, Float64, order; space=:P)
M = FESpace(Γ, reffe_M; conformity=:L2, dirichlet_tags="boundary")
L = FESpaces.PolytopalFESpace(Ω, Float64, order+1; space=:P)
N = TrialFESpace(M,g)

mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
X = MultiFieldFESpace([V, N];style=mfs)
Y = MultiFieldFESpace([V, M];style=mfs)

x = randn(num_free_dofs(X))
y = randn(num_free_dofs(Y))
xh = FEFunction(X,x)
yh = FEFunction(Y,y) # Zero Dirichlet values

# Patch operator, multi-field

R = reconstruction_operator(ptopo,L,Y,Ωp,Γp,dΩp,dΓp)
@test isa(R,LocalOperator{<:Geometry.PatchTopology})

RX = R(X)
RY = R(Y)
@test isa(RX,MultiFieldFESpace)
@test MultiFieldStyle(RX) == MultiFieldStyle(X)
@test num_free_dofs(RX) == num_free_dofs(X)
@test all(map((a,b) -> num_dirichlet_dofs(a) == num_dirichlet_dofs(b), RX, X))
@test all(sf -> get_triangulation(sf) === Ω, RX) # trian_out of R

# The dofs are those of PatchFESpace
Xp = FESpaces.PatchFESpace(X,ptopo)
@test all(map((a,b) -> get_cell_dof_ids(a) == get_cell_dof_ids(b), RX, Xp))

# Functions: FEFunction(R(X),x) == R(xh), including the Dirichlet data
RxT, RxF = FEFunction(RX,x)
Rxh = R(xh)
@test l2(Rxh,dΩp) > 0.1
@test l2(RxT + RxF - Rxh,dΩp) < 1.e-10

# Bases: the bases of R(X) are R applied to the bases of X
a((uT,uF),(vT,vF)) = ∫(∇(uT+uF)⋅∇(vT+vF))dΩp
assem = FESpaces.PatchAssembler(ptopo,X,Y)
mats_old = assemble_matrix(assem,FESpaces.collect_patch_cell_matrix(
  assem,Xp,Xp,a(R(get_trial_fe_basis(X)),R(get_fe_basis(Y)))
))
mats_new = assemble_matrix(assem,FESpaces.collect_patch_cell_matrix(
  assem,Xp,Xp,a(get_trial_fe_basis(RX),get_fe_basis(RY))
))
@test all(p -> all(mats_old[p].array .≈ mats_new[p].array), 1:Geometry.num_patches(ptopo))

# Global assembly on R(X), with a standard assembler
op = AffineFEOperator(a,v -> DomainContribution(),RX,RY)
A, b = get_matrix(op), get_vector(op)
Ryh = R(yh)
@test dot(y,A*x-b) ≈ sum(∫(∇(Rxh)⋅∇(Ryh))dΩp)

# The output space as input of another (cell-wise) operator
PΩ_mf = mf_projection_operator(V,Ωp,dΩp)
PRX = PΩ_mf(RX)
PRxT, PRxF = FEFunction(PRX,x)
@test l2(PRxT + PRxF - PΩ_mf(FEFunction(RX,x)),dΩp) < 1.e-10

# Cell-wise operator, single-field

xΓ = get_free_dof_values(xh[2])
PΓ = projection_operator(M,Γp,dΓp)
@test isa(PΓ,LocalOperator{Nothing})

PN = PΓ(N)
@test isa(PN,SingleFieldFESpace)
@test get_triangulation(PN) === Γp
@test get_cell_dof_ids(PN) == get_cell_dof_ids(N,Γp)
@test l2(FEFunction(PN,xΓ) - xh[2],dΓp) < 1.e-10 # P is the identity on M

# Cell-wise operator, multi-field

PΓ_mf = mf_projection_operator(M,Γp,dΓp)
PX = PΓ_mf(X)
@test isa(PX,MultiFieldFESpace)
PxT, PxF = FEFunction(PX,x)
@test l2(PxT + PxF - PΓ_mf(xh),dΓp) < 1.e-10

end # module
