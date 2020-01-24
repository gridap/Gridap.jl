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
@test isa(fun_Γ,SkeletonPair)

end # module
