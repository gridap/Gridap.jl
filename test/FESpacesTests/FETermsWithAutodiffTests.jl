module FETermsWithAutodiff

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

res(uh,v) = ∇(uh)⋅∇(v)
jac(uh,u,v) = ∇(u)⋅∇(v)

t_Ω = FETerm(res,jac,trian,quad)
t_auto_Ω = FETerm(res,trian,quad)

cell_r = get_cell_residual(t_Ω,uh,dv)
cell_j = get_cell_jacobian(t_Ω,uh,du,dv)

cell_r_auto = get_cell_residual(t_auto_Ω,uh,dv)
cell_j_auto = get_cell_jacobian(t_auto_Ω,uh,du,dv)

test_array(cell_r_auto,cell_r)
test_array(cell_j_auto,cell_j)

trian_Γ = BoundaryTriangulation(model)
quad_Γ = CellQuadrature(trian_Γ,2)

t_Γ = FETerm(res,jac,trian_Γ,quad_Γ)
t_auto_Γ = FETerm(res,trian_Γ,quad_Γ)

cell_r_Γ = get_cell_residual(t_Γ,uh,dv)
cell_j_Γ = get_cell_jacobian(t_Γ,uh,du,dv)

cell_r_Γ_auto = get_cell_residual(t_auto_Γ,uh,dv)
cell_j_Γ_auto = get_cell_jacobian(t_auto_Γ,uh,du,dv)

test_array(cell_r_Γ_auto,cell_r_Γ)
test_array(cell_j_Γ_auto,cell_j_Γ)

end # module
