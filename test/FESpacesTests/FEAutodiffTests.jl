module FEAutodiffTests

using LinearAlgebra
using Gridap.FESpaces
using Gridap.Arrays
using Gridap.Fields
using Gridap.Geometry

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

V = FESpace(model=model,valuetype=Float64,reffe=:Lagrangian,order=1,conformity=:H1)
U = TrialFESpace(V)

dv = get_cell_basis(V)
du = get_cell_basis(U)
uh = FEFunction(U,rand(num_free_dofs(U)))

trian = Triangulation(model)
quad = CellQuadrature(trian,2)

function user_uh_to_cell_energy(uh)
  cell_e = integrate(0.5*∇(uh)⋅∇(uh),trian,quad)
end

function user_uh_to_cell_residual(uh)
  cell_r = integrate(∇(uh)⋅∇(dv),trian,quad)
end

function user_uh_to_cell_jacobian(uh)
  cell_j = integrate(∇(du)⋅∇(dv),trian,quad)
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

end # module
