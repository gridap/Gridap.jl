module FEAutodiffTests

using LinearAlgebra
using Gridap.FESpaces
using Gridap.Arrays
using Gridap.Fields
using Gridap.Geometry
using Gridap.TensorValues 
using Gridap.CellData
using Gridap.ReferenceFEs

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

Ω = Triangulation(model)
dΩ = LebesgueMeasure(Ω,2)

V = FESpace(model,ReferenceFE(:Lagrangian,Float64,2),conformity=:H1)
U = V

dv = get_cell_shapefuns(V)
du = get_cell_shapefuns_trial(U)
uh = FEFunction(U,rand(num_free_dofs(U)))

function user_uh_to_cell_energy(uh)
  cell_e = ∫( 0.5*∇(uh)⋅∇(uh) )*dΩ
end

function user_uh_to_cell_residual(uh)
  cell_r = ∫(∇(uh)⋅∇(dv))*dΩ
end

function user_uh_to_cell_jacobian(uh)
  cell_j = ∫(∇(du)⋅∇(dv))*dΩ
end

cell_r = get_array(user_uh_to_cell_residual(uh))
cell_j = get_array(user_uh_to_cell_jacobian(uh))
cell_h = cell_j

cell_r_auto = autodiff_cell_residual_from_energy(user_uh_to_cell_energy,uh)
cell_h_auto = autodiff_cell_jacobian_from_energy(user_uh_to_cell_energy,uh)
cell_j_auto = autodiff_cell_jacobian_from_residual(user_uh_to_cell_residual,uh)

test_array(cell_r_auto,cell_r)
test_array(cell_j_auto,cell_j)
test_array(cell_h_auto,cell_h)

Γ = BoundaryTriangulation(model)
dΓ = LebesgueMeasure(Γ,2)
cell_ids = get_cell_id(Γ)

function user_uh_to_cell_energy_Γ(uh)
  cell_e = ∫( 0.5*∇(uh)⋅∇(uh) )*dΓ
end

function user_uh_to_cell_residual_Γ(uh)
  cell_r = ∫( ∇(uh)⋅∇(dv) )*dΓ
end

function user_uh_to_cell_jacobian_Γ(uh)
  cell_j = ∫( ∇(du)⋅∇(dv) )*dΓ
end

cell_e_Γ = get_array(user_uh_to_cell_energy_Γ(uh))
cell_r_Γ = get_array(user_uh_to_cell_residual_Γ(uh))
cell_j_Γ = get_array(user_uh_to_cell_jacobian_Γ(uh))
cell_h_Γ = cell_j_Γ

cell_r_Γ_auto = autodiff_cell_residual_from_energy(user_uh_to_cell_energy_Γ,uh,cell_ids)
cell_h_Γ_auto = autodiff_cell_jacobian_from_energy(user_uh_to_cell_energy_Γ,uh,cell_ids)
cell_j_Γ_auto = autodiff_cell_jacobian_from_residual(user_uh_to_cell_residual_Γ,uh,cell_ids)

test_array(cell_r_Γ_auto,cell_r_Γ)
test_array(cell_j_Γ_auto,cell_j_Γ)
test_array(cell_h_Γ_auto,cell_h_Γ)

const p = 3
j(∇u) = norm(∇u)^(p-2) * ∇u
dj(∇du,∇u) = (p-2)*norm(∇u)^(p-4)*inner(∇u,∇du)*∇u + norm(∇u)^(p-2)*∇du
f(x) = 0

res(u,v) = ∫( ∇(v)⋅(j∘∇(u)) - v*f)*dΩ
jac(u,du,v) = ∫( ∇(v)⋅(dj∘(∇(du),∇(u))) )*dΩ

function user_uh_to_cell_residual_2(uh)
  cell_r = res(uh,dv)
end

function user_uh_to_cell_jacobian_2(uh)
  cell_j = jac(uh,du,dv)
end

cell_j = get_array(user_uh_to_cell_jacobian_2(uh))

cell_j_auto = autodiff_cell_jacobian_from_residual(user_uh_to_cell_residual_2,uh)

test_array(cell_j_auto,cell_j)

end # module
