using Gridap
using Gridap.Algebra
using Gridap.Arrays
using Gridap.Geometry
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.MultiField

model = CartesianDiscreteModel((0,1,0,1),(3,3))
Ω1 = Triangulation(model,[1,2,3])
Λ = SkeletonTriangulation(model)

V = FESpace(Ω1,ReferenceFE(lagrangian,Float64,1))
uh = zero(V);
dΛ = Measure(Λ,2)

f2(xh,yh) = ∫(mean(xh)*mean(yh))dΛ
dv = get_fe_basis(V);
j = jacobian(uh->f2(uh,dv),uh);
J = assemble_matrix(j,V,V)