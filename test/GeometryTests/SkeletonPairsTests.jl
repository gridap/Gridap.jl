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

a = rand(num_cells(trian))
a_Γ = reindex(a,trian_Γ)
@test isa(a_Γ,SkeletonPair)


end # module
