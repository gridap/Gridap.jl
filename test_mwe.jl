using Gridap
using Gridap.Algebra
using Gridap.Arrays
using Gridap.Geometry
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.MultiField

using Gridap.Arrays: lazy_collect
using FillArrays

model = CartesianDiscreteModel((0,1,0,1),(3,3))
Ω1 = Triangulation(model,[1,2,3])
Λ = SkeletonTriangulation(model)

V = FESpace(Ω1,ReferenceFE(lagrangian,Float64,1))
uh = zero(V);
dΛ = Measure(Λ,2)

f2(xh,yh) = ∫(mean(xh)*mean(yh))dΛ
dv = get_fe_basis(V);
j = jacobian(uh->f2(uh,dv),uh)
J = assemble_matrix(j,V,V)

bmats = j[Λ]
lazy_collect(bmats)

a_plus = bmats.args[1]
a_minus = bmats.args[2]
lazy_collect(a_plus)
lazy_collect(a_minus)


j_to_cfg, j_to_ydual = a_plus.args
lazy_collect(j_to_cfg)
lazy_collect(j_to_ydual)

n = length(j_to_cfg)
j_to_result = map(return_cache,Fill(AutoDiffMap(),n),j_to_cfg,j_to_ydual)


