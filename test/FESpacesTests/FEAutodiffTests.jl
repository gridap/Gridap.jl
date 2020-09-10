module FEAutodiffTests

using LinearAlgebra
using Gridap.FESpaces
using Gridap.Arrays
using Gridap.Fields
using Gridap.Geometry
using Gridap.TensorValues 
using Gridap.CellData

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

V = FESpace(model=model,valuetype=Float64,reffe=:Lagrangian,order=2,conformity=:H1)
U = TrialFESpace(V)

dv = get_cell_basis(V)
du = get_cell_basis(U)
uh = FEFunction(U,rand(num_free_dofs(U)))

trian = Triangulation(model)
quad = CellQuadrature(trian,2)

function user_uh_to_cell_energy(uh)
  cell_e = integrate(0.5*∇(uh)⋅∇(uh),quad)
end

function user_uh_to_cell_residual(uh)
  cell_r = integrate(∇(uh)⋅∇(dv),quad)
end

function user_uh_to_cell_jacobian(uh)
  cell_j = integrate(∇(du)⋅∇(dv),quad)
end

cell_r = user_uh_to_cell_residual(uh)
cell_j = user_uh_to_cell_jacobian(uh)
cell_h = cell_j

cell_r_auto = autodiff_cell_residual_from_energy(user_uh_to_cell_energy,uh)
cell_h_auto = autodiff_cell_jacobian_from_energy(user_uh_to_cell_energy,uh)
cell_j_auto = autodiff_cell_jacobian_from_residual(user_uh_to_cell_residual,uh)

test_array(cell_r_auto,cell_r)
test_array(cell_j_auto,cell_j)
test_array(cell_h_auto,cell_h)

trian_Γ = BoundaryTriangulation(model)
quad_Γ = CellQuadrature(trian_Γ,2)
cell_ids = get_cell_id(trian_Γ)

function user_uh_to_cell_energy_Γ(uh)
  cell_e = integrate(0.5*∇(uh)⋅∇(uh),quad_Γ)
end

function user_uh_to_cell_residual_Γ(uh)
  cell_r = integrate(∇(uh)⋅∇(dv),quad_Γ)
end

function user_uh_to_cell_jacobian_Γ(uh)
  cell_j = integrate(∇(du)⋅∇(dv),quad_Γ)
end

cell_e_Γ = user_uh_to_cell_energy_Γ(uh)
cell_r_Γ = user_uh_to_cell_residual_Γ(uh)
cell_j_Γ = user_uh_to_cell_jacobian_Γ(uh)
cell_h_Γ = cell_j_Γ

cell_r_Γ_auto = autodiff_cell_residual_from_energy(user_uh_to_cell_energy_Γ,uh,cell_ids)
cell_h_Γ_auto = autodiff_cell_jacobian_from_energy(user_uh_to_cell_energy_Γ,uh,cell_ids)
cell_j_Γ_auto = autodiff_cell_jacobian_from_residual(user_uh_to_cell_residual_Γ,uh,cell_ids)

test_array(cell_r_Γ_auto,cell_r_Γ)
test_array(cell_j_Γ_auto,cell_j_Γ)
test_array(cell_h_Γ_auto,cell_h_Γ)

const p = 3
@law j(∇u) = norm(∇u)^(p-2) * ∇u
@law dj(∇du,∇u) = (p-2)*norm(∇u)^(p-4)*inner(∇u,∇du)*∇u + norm(∇u)^(p-2)*∇du
f(x) = 0

res(u,v) = ∇(v)⋅j(∇(u)) - v*f
jac(u,du,v) = ∇(v)⋅dj(∇(du),∇(u))

function user_uh_to_cell_residual_2(uh)
  cell_r = integrate(res(uh,dv),quad)
end

function user_uh_to_cell_jacobian_2(uh)
  cell_j = integrate(jac(uh,du,dv),quad)
end

cell_j = user_uh_to_cell_jacobian_2(uh)

cell_j_auto = autodiff_cell_jacobian_from_residual(user_uh_to_cell_residual_2,uh)

test_array(cell_j_auto,cell_j)

end # module
