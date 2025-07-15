using Gridap, Gridap.FESpaces, Gridap.MultiField, Gridap.CellData, Gridap.Helpers, Gridap.Fields
using Gridap.FESpaces: _change_argument, _compute_cell_ids, autodiff_array_gradient
using Gridap.MultiField: blocks, MultiFieldFEFunction, mortar
using Gridap.CellData: is_change_possible
using Test

domain = (0,1,0,1)
partition = (10,10)
model = CartesianDiscreteModel(domain,partition)
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
dΩ = Measure(Ω,2)
dΓ = Measure(Γ,2)
V1 = FESpace(Γ,ReferenceFE(lagrangian,Float64,1))
# V1 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
V2 = FESpace(model,ReferenceFE(lagrangian,Float64,2))
X = MultiFieldFESpace([V1,V2])
f(xh) = ∫(xh[1]+xh[2])dΓ
uh = zero(X)
du = gradient(f,uh)
du1 = gradient(x->f((x,uh[2])),uh[1])
du2 = gradient(x->f((uh[1],x)),uh[2])

@test lazy_map(Gridap.MultiField.GetIndex(1),du[Γ]) == du1[Γ]
@test lazy_map(Gridap.MultiField.GetIndex(2),du[Γ]) == du2[Γ]

du1_vec = assemble_vector(du1,V1)
du2_vec = assemble_vector(du2,V2)
du_vec = assemble_vector(du,X)

@test du_vec == [du1_vec;du2_vec]

f2(xh,yh) = ∫(xh[1]⋅yh[1]+xh[2]⋅yh[2]+xh[1]⋅xh[2]⋅yh[1])dΓ
dv = get_fe_basis(X)
j = jacobian(uh->f2(uh,dv),uh)
J = assemble_matrix(j,X,X)

f2_jac(xh,dxh,yh) = ∫(dxh[1]⋅yh[1]+dxh[2]⋅yh[2]+dxh[1]⋅xh[2]⋅yh[1]+xh[1]⋅dxh[2]⋅yh[1])dΓ
op = FEOperator(f2,f2_jac,X,X)
J_fwd = jacobian(op,uh)

@test J_fwd == J

## Hessian

e((uh,ph)) = ∫( uh*uh + uh*ph + ph*ph )dΓ
# dv = get_fe_basis(X)
j_hes = hessian(e,uh)
J_hes = assemble_matrix(j_hes,X,X)
