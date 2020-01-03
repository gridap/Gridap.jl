module SkeletonPairsTests

using Test
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

trian_Γ = SkeletonTriangulation(model)
trian = get_volume_triangulation(trian_Γ)

fun(x) = sin(pi*x[1])*cos(pi*x[2])
q2x = get_cell_map(trian)
funq = compose(fun,q2x)
fun_Γ = restrict(funq,trian_Γ)

s = CompressedArray([Point{1,Float64}[(0.25,),(0.75,)]],get_cell_type(trian_Γ))
j = jump(fun_Γ)
jx = evaluate(j,s)
r = fill([0.0,0.0],size(jx))
test_array(jx,r)

shapefuns = get_cell_shapefuns(trian)
shapefuns_Γ = restrict(shapefuns,trian_Γ)
j = jump(shapefuns_Γ)
jx = evaluate(j,s)
r = collect(jx.left)
test_array(jx.left,r)
r = collect(jx.right)
test_array(jx.right,r)

jac = ∇(get_cell_map(trian))
jac_Γ = restrict(jac,trian_Γ)
j = jump(jac_Γ)
jx = evaluate(j,s)
r = collect(jx)
test_array(jx,r)

end # module
