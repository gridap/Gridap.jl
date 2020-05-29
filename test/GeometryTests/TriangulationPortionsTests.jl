module TriangulationPortionsTests

using Gridap.ReferenceFEs
using Gridap.Geometry
using Test

domain = (0,1,0,1)
partition = (10,10)
oldmodel = CartesianDiscreteModel(domain,partition)
oldtrian = get_triangulation(oldmodel)

cell_to_oldcell = collect(1:34)
trian = TriangulationPortion(oldtrian,cell_to_oldcell)
test_triangulation(trian)

oldbtrian = BoundaryTriangulation(oldmodel)

bface_to_oldbface = collect(1:32)
btrian = TriangulationPortion(oldbtrian,bface_to_oldbface)
test_triangulation(btrian)

nb = get_normal_vector(btrian)
@test length(nb) == num_cells(btrian)

sface_to_oldsface = collect(10:45)
oldstrian = SkeletonTriangulation(oldmodel)
strian = TriangulationPortion(oldstrian,sface_to_oldsface)
test_triangulation(strian)

ns = get_normal_vector(strian)
@test length(ns) == num_cells(strian)

v = [i for i=1:num_cells(trian)]
cell_to_oldcell = [2,9,7]
trian_portion = TriangulationPortion(oldtrian,cell_to_oldcell)
trian_portion_portion = TriangulationPortion(trian_portion,[3,1])
@test get_cell_id(trian_portion_portion)==[7,2]
@test restrict(v, trian_portion_portion)==[7,2]

#using Gridap.Visualization
#writevtk(oldtrian,"oldtrian")
#writevtk(btrian,"btrian",cellfields=["normal"=>nb],celldata=["oldcell"=>get_cell_id(btrian)])
#writevtk(strian,"strian",cellfields=["normal"=>ns],
#  celldata=["oldcell_left"=>get_cell_id(strian).left,"oldcell_right"=>get_cell_id(strian).right])

end # module
