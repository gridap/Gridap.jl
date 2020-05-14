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

#using Gridap.Visualization
#writevtk(oldtrian,"oldtrian")
#writevtk(btrian,"btrian",cellfields=["normal"=>nb],celldata=["oldcell"=>get_cell_id(btrian)])

end # module

