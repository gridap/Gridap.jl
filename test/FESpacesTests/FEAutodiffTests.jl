module FEAutodiffTests

using Test
using LinearAlgebra
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.Arrays
using Gridap.Fields
using Gridap.Geometry
using Gridap.TensorValues
using Gridap.CellData
using Gridap.ReferenceFEs
using ForwardDiff
using SparseArrays

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

Ω = Triangulation(model)
dΩ = Measure(Ω,2)

V = FESpace(model,ReferenceFE(lagrangian,Float64,2),conformity=:H1)
U = V

dv = get_fe_basis(V)
du = get_trial_fe_basis(U)
uh = FEFunction(U,rand(num_free_dofs(U)))

ener(uh) = ∫( 0.5*∇(uh)⋅∇(uh) )*dΩ
res(uh) = ∫(∇(uh)⋅∇(dv))*dΩ
jac(uh) = ∫(∇(du)⋅∇(dv))*dΩ

cell_r = get_array(res(uh))
cell_j = get_array(jac(uh))
cell_h = cell_j

cell_r_auto = get_array(gradient(ener,uh))
cell_j_auto = get_array(jacobian(res,uh))
cell_h_auto = get_array(hessian(ener,uh))

test_array(cell_r_auto,cell_r,≈)
test_array(cell_j_auto,cell_j,≈)
test_array(cell_h_auto,cell_h,≈)

Γ = BoundaryTriangulation(model)
dΓ = Measure(Γ,2)

ener(uh) = ∫( 0.5*∇(uh)⋅∇(uh) )*dΓ
res(uh) = ∫( ∇(uh)⋅∇(dv) )*dΓ
jac(uh) = ∫( ∇(du)⋅∇(dv) )*dΓ

cell_r = get_array(res(uh))
cell_j = get_array(jac(uh))
cell_h = cell_j

cell_r_auto = get_array(gradient(ener,uh))
cell_j_auto = get_array(jacobian(res,uh))
cell_h_auto = get_array(hessian(ener,uh))

test_array(cell_r_auto,cell_r,≈)
test_array(cell_j_auto,cell_j,≈)
test_array(cell_h_auto,cell_h,≈)

ener(uh) = ∫( 0.5*∇(uh)⋅∇(uh) )*dΓ + ∫( 0.5*∇(uh)⋅∇(uh) )*dΩ
res(uh) = ∫( ∇(uh)⋅∇(dv) )*dΓ + ∫(∇(uh)⋅∇(dv))*dΩ
jac(uh) = ∫( ∇(du)⋅∇(dv) )*dΓ + ∫(∇(du)⋅∇(dv))*dΩ

cell_r = res(uh)
cell_j = jac(uh)
cell_h = cell_j

cell_r_auto = gradient(ener,uh)
cell_j_auto = jacobian(res,uh)
cell_h_auto = hessian(ener,uh)

test_array(cell_r_auto[Ω],cell_r[Ω],≈)
test_array(cell_j_auto[Ω],cell_j[Ω],≈)
test_array(cell_h_auto[Ω],cell_h[Ω],≈)

test_array(cell_r_auto[Γ],cell_r[Γ],≈)
test_array(cell_j_auto[Γ],cell_j[Γ],≈)
test_array(cell_h_auto[Γ],cell_h[Γ],≈)

const p = 3
j(∇u) = norm(∇u)^(p-2) * ∇u
dj(∇du,∇u) = (p-2)*norm(∇u)^(p-4)*inner(∇u,∇du)*∇u + norm(∇u)^(p-2)*∇du
f(x) = 0

res(u,v) = ∫( ∇(v)⋅(j∘∇(u)) - v*f)*dΩ
jac(u,du,v) = ∫( ∇(v)⋅(dj∘(∇(du),∇(u))) )*dΩ

cell_j = get_array(jac(uh,du,dv))
cell_j_auto = get_array(jacobian(u->res(u,dv),uh))

test_array(cell_j_auto,cell_j,≈)

## comparing AD of integration over Skeleton faces with ForwardDiff results ##

model = CartesianDiscreteModel((0.,1.,0.,1.),(2,2))
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
Λ = SkeletonTriangulation(model)

dΩ = Measure(Ω,2)
dΓ = Measure(Γ,2)
dΛ = Measure(Λ,2)

n_Γ = get_normal_vector(Γ)
n_Λ = get_normal_vector(Λ)

reffe = ReferenceFE(lagrangian,Float64,1)
V = TestFESpace(model,reffe,conformity=:L2)

u(x) = sin(norm(x))
U = TrialFESpace(V)

uh = FEFunction(U,rand(num_free_dofs(U)))

f_Λ(uh) = ∫(mean(uh))*dΛ
a_Λ(u) = ∫( - jump(u*n_Λ)⊙mean(∇(u))
            - mean(Δ(u))
            + jump(u*n_Λ)⊙jump(u*n_Λ) )dΛ

# functionals having mean and jump of products of FEFunctions/CellFields
g_Λ(uh) = ∫(mean(uh*uh))*dΛ
h_Λ(uh) = ∫(jump(uh*uh*n_Λ)⋅jump(uh*uh*n_Λ))*dΛ
j_Λ(uh) = ∫(mean(∇(uh)*uh)⋅jump((∇(uh)⋅∇(uh))*n_Λ))*dΛ

@test sum(∫(mean(uh*uh))*dΛ) == 0.5*sum(∫(uh.plus*uh.plus + uh.minus*uh.minus)*dΛ)

function f_uh_free_dofs(f,uh,θ)
  dir = similar(uh.dirichlet_values,eltype(θ))
  uh = FEFunction(U,θ,dir)
  sum(f(uh))
end

f_Λ_(θ) = f_uh_free_dofs(f_Λ,uh,θ)
a_Λ_(θ) = f_uh_free_dofs(a_Λ,uh,θ)
θ = get_free_dof_values(uh)

# testing gradients of functionals involving Skeleton integration terms

gridapgradf = assemble_vector(gradient(f_Λ,uh),U)
fdgradf = ForwardDiff.gradient(f_Λ_,θ)
test_array(gridapgradf,fdgradf,≈)

gridapgrada = assemble_vector(gradient(a_Λ,uh),U)
fdgrada = ForwardDiff.gradient(a_Λ_,θ)
test_array(gridapgrada,fdgrada,≈)

g_Λ_(θ) = f_uh_free_dofs(g_Λ,uh,θ)
h_Λ_(θ) = f_uh_free_dofs(h_Λ,uh,θ)
j_Λ_(θ) = f_uh_free_dofs(j_Λ,uh,θ)

gridapgradg = assemble_vector(gradient(g_Λ,uh),U)
fdgradg = ForwardDiff.gradient(g_Λ_,θ)
test_array(gridapgradg,fdgradg,≈)

gridapgradh = assemble_vector(gradient(h_Λ,uh),U)
fdgradh = ForwardDiff.gradient(h_Λ_,θ)
test_array(gridapgradh,fdgradh,≈)

gridapgradj = assemble_vector(gradient(j_Λ,uh),U)
fdgradj = ForwardDiff.gradient(j_Λ_,θ)
test_array(gridapgradj,fdgradj,≈)

# testing jacobians of functionals involving Skeleton inetgration terms
uh = FEFunction(U,rand(num_free_dofs(U)))
θ = get_free_dof_values(uh)
dv = get_fe_basis(V)
du = get_trial_fe_basis(U)

a(uh,vh) = ∫(mean(uh)*mean(vh))*dΛ
# b(uh,vh) = ∫(-jump(uh*n_Λ)⊙mean(∇(vh)) +
#               jump(vh*n_Λ)⊙mean(∇(uh)) +
#               mean(Δ(uh))*mean(Δ(vh))  +
#               jump(uh*n_Λ)⊙jump(vh*n_Λ))*dΛ
# g(uh,vh) = ∫(mean(uh*vh))*dΛ
h(uh,vh) = ∫(jump(uh*vh*n_Λ)⊙jump(vh*vh*n_Λ))*dΛ
# j(uh,vh) = ∫(mean(∇(vh)*uh)⊙jump((∇(vh)⊙∇(uh))*n_Λ))*dΛ

# We use the strongly-typed lowel-level interface of SparseMatrixCSC{T,Int} here
# as ForwardDiff doesn't work directly through `assemble_vector``, this is due
# to the typing of SparseMatrixCSC{T,Int} and Vector{T} inside high-level
# `assemble_vector` API (via the `SparseMatrixAssembler`) using the type T of
# the eltype of FESpace U, which remains Float even though Dual numbers are
# passed into θ.
# So as to mitigate this problem, low-level interface of assemble_vector is
# being used with SparseMatrixCSC{T,Int} and Vector{T} constructed by hand with
# types of θ which can be dual numbers when using ForwardDiff or even #
# ReverseDiff types, when using ReverseDiff! This is just for testing purposes
function _change_input(f,θ,uh)
  dir = similar(uh.dirichlet_values,eltype(θ))
  U = uh.fe_space
  uh = FEFunction(U,θ,dir)
  T = eltype(θ)
  matrix_type = SparseMatrixCSC{T,Int}
  vector_type = Vector{T}
  assem = SparseMatrixAssembler(matrix_type,vector_type,U,U)
  assemble_vector(f(uh),assem,U)
end

function _assemble_jacobian(f,uh)
  assemble_matrix(jacobian(f,uh),U,V)
end

f(uh) = a(uh,dv)
jac_gridap_a = _assemble_jacobian(f,uh)
collect(get_array(jacobian(f,uh))) # just to check the working
f_(θ) = _change_input(f,θ,uh)
jac_forwdiff_a = ForwardDiff.jacobian(f_,θ)
test_array(jac_gridap_a,jac_forwdiff_a,≈)

# f(uh) = b(uh,dv)
# jac_gridap_b = _assemble_jacobian(f,uh)
# f_(θ) = _change_input(f,θ,uh)
# jac_forwdiff_b = ForwardDiff.jacobian(f_,θ)
# test_array(jac_gridap_b,jac_forwdiff_b,≈)

# f(uh) = g(uh,dv)
# jac_gridap_g = _assemble_jacobian(f,uh)
# f_(θ) = _change_input(f,θ,uh)
# jac_forwdiff_g = ForwardDiff.jacobian(f_,θ)
# test_array(jac_gridap_g,jac_forwdiff_g,≈)

f(uh) = h(uh,dv)
jac_gridap_h = _assemble_jacobian(f,uh)
f_(θ) = _change_input(f,θ,uh)
jac_forwdiff_h = ForwardDiff.jacobian(f_,θ)
test_array(jac_gridap_h,jac_forwdiff_h,≈)

# f(uh) = j(uh,dv)
# jac_gridap_j = _assemble_jacobian(f,uh)
# f_(θ) = _change_input(f,θ,uh)
# jac_forwdiff_j = ForwardDiff.jacobian(f_,θ)
# test_array(jac_gridap_j,jac_forwdiff_j,≈)

end # module
