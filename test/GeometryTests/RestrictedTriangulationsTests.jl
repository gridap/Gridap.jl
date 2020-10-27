module RestrictedTriangulationsTests

using Gridap.ReferenceFEs
using Gridap.Geometry
using Test

domain = (0,1,0,1)
partition = (10,10)
oldmodel = CartesianDiscreteModel(domain,partition)
oldtrian = get_triangulation(oldmodel)

cell_to_oldcell = collect(1:34)
trian = RestrictedTriangulation(oldtrian,cell_to_oldcell)
test_triangulation(trian)

cell_to_oldcell = [2,9,7]
trian_portion = RestrictedTriangulation(oldtrian,cell_to_oldcell)
trian_portion_portion = RestrictedTriangulation(trian_portion,[3,1])
test_triangulation(trian_portion_portion)
@test get_cell_id(trian_portion_portion)==[7,2]

#using Gridap.Visualization
#writevtk(oldtrian,"oldtrian")
#writevtk(btrian,"btrian",cellfields=["normal"=>nb],celldata=["oldcell"=>get_cell_id(btrian)])
#writevtk(strian,"strian",cellfields=["normal"=>ns],
#  celldata=["oldcell_left"=>get_cell_id(strian).left,"oldcell_right"=>get_cell_id(strian).right])

end # module
