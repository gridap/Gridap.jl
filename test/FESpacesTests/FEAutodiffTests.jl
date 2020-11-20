module FEAutodiffTests

using LinearAlgebra
using Gridap.Algebra
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

ener(uh) = ∫( 0.5*∇(uh)⋅∇(uh) )*dΩ
res(uh) = ∫(∇(uh)⋅∇(dv))*dΩ
jac(uh) = ∫(∇(du)⋅∇(dv))*dΩ

cell_r = get_array(res(uh))
cell_j = get_array(jac(uh))
cell_h = cell_j

cell_r_auto = get_array(gradient(ener,uh))
cell_j_auto = get_array(jacobian(res,uh))
cell_h_auto = get_array(hessian(ener,uh))

test_array(cell_r_auto,cell_r)
test_array(cell_j_auto,cell_j)
test_array(cell_h_auto,cell_h)

Γ = BoundaryTriangulation(model)
dΓ = LebesgueMeasure(Γ,2)

ener(uh) = ∫( 0.5*∇(uh)⋅∇(uh) )*dΓ
res(uh) = ∫( ∇(uh)⋅∇(dv) )*dΓ
jac(uh) = ∫( ∇(du)⋅∇(dv) )*dΓ

cell_r = get_array(res(uh))
cell_j = get_array(jac(uh))
cell_h = cell_j

cell_r_auto = get_array(gradient(ener,uh))
cell_j_auto = get_array(jacobian(res,uh))
cell_h_auto = get_array(hessian(ener,uh))

test_array(cell_r_auto,cell_r)
test_array(cell_j_auto,cell_j)
test_array(cell_h_auto,cell_h)

ener(uh) = ∫( 0.5*∇(uh)⋅∇(uh) )*dΓ + ∫( 0.5*∇(uh)⋅∇(uh) )*dΩ
res(uh) = ∫( ∇(uh)⋅∇(dv) )*dΓ + ∫(∇(uh)⋅∇(dv))*dΩ
jac(uh) = ∫( ∇(du)⋅∇(dv) )*dΓ + ∫(∇(du)⋅∇(dv))*dΩ

cell_r = res(uh)
cell_j = jac(uh)
cell_h = cell_j

cell_r_auto = gradient(ener,uh)
cell_j_auto = jacobian(res,uh)
cell_h_auto = hessian(ener,uh)

test_array(cell_r_auto[Ω],cell_r[Ω])
test_array(cell_j_auto[Ω],cell_j[Ω])
test_array(cell_h_auto[Ω],cell_h[Ω])

test_array(cell_r_auto[Γ],cell_r[Γ])
test_array(cell_j_auto[Γ],cell_j[Γ])
test_array(cell_h_auto[Γ],cell_h[Γ])

const p = 3
j(∇u) = norm(∇u)^(p-2) * ∇u
dj(∇du,∇u) = (p-2)*norm(∇u)^(p-4)*inner(∇u,∇du)*∇u + norm(∇u)^(p-2)*∇du
f(x) = 0

res(u,v) = ∫( ∇(v)⋅(j∘∇(u)) - v*f)*dΩ
jac(u,du,v) = ∫( ∇(v)⋅(dj∘(∇(du),∇(u))) )*dΩ

cell_j = get_array(jac(uh,du,dv))
cell_j_auto = get_array(jacobian(u->res(u,dv),uh))

test_array(cell_j_auto,cell_j)

end # module
