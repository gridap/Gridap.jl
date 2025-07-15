using Gridap, Gridap.FESpaces, Gridap.MultiField, Gridap.CellData, Gridap.Helpers, Gridap.Fields
using Gridap.FESpaces: _change_argument, _compute_cell_ids, autodiff_array_gradient
using Gridap.MultiField: blocks, MultiFieldFEFunction, mortar
using Gridap.CellData: is_change_possible
using Test

domain = (0,1,0,1)
partition = (20,20)
model = CartesianDiscreteModel(domain,partition)
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
dΩ = Measure(Ω,2)
dΓ = Measure(Γ,2)
# V1 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
# V2 = FESpace(Γ,ReferenceFE(lagrangian,Float64,2))
# X = MultiFieldFESpace([V1,V2])
# f(xh) = ∫(xh[1]+xh[2])dΓ
# uh = zero(X)
# du = gradient(f,uh)
# du1 = gradient(x->f((x,uh[2])),uh[1])
# du2 = gradient(x->f((uh[1],x)),uh[2])

# @test lazy_map(Gridap.MultiField.GetIndex(1),du[Γ]) == du1[Γ]
# @test lazy_map(Gridap.MultiField.GetIndex(2),du[Γ]) == du2[Γ]

# du1_vec = assemble_vector(du1,V1)
# du2_vec = assemble_vector(du2,V2)
# du_vec = assemble_vector(du,X)

# @test du_vec == [du1_vec;du2_vec]

# Timing
using BenchmarkTools
V1 = FESpace(Ω,ReferenceFE(lagrangian,Float64,1))
V2 = FESpace(Ω,ReferenceFE(lagrangian,Float64,2))
V3 = FESpace(Ω,ReferenceFE(lagrangian,Float64,3))
V4 = FESpace(Ω,ReferenceFE(lagrangian,Float64,4))
X = MultiFieldFESpace([V1,V2,V3,V4])
f(xh) = ∫(xh[1]+xh[2]+xh[3]+xh[4])dΓ
uh = zero(X)
@benchmark gradient($f,$uh)
grad = gradient(f,uh)
@benchmark assemble_vector($grad,$X)

@benchmark MultiField.my_new_gradient($f,$uh)
grad = MultiField.my_new_gradient(f,uh)
@benchmark assemble_vector($grad,$X)