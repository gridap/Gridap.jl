module RestrictedTriangulationsTests

using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Test

domain = (-1,1,-1,1)
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

const R = 0.7

function is_in(coords)
  n = length(coords)
  x = (1/n)*sum(coords)
  d = x[1]^2 + x[2]^2 - R^2
  d < 0
end

oldcell_to_coods = get_cell_coordinates(oldtrian)
oldcell_to_mask = lazy_map(is_in,oldcell_to_coods)
trian = RestrictedTriangulation(oldtrian,oldcell_to_mask)
trian = Triangulation(oldmodel,oldcell_to_mask)
trian = Triangulation(oldtrian,oldcell_to_mask)
trian = Triangulation(oldmodel,tags="interior")
trian = Triangulation(oldmodel,get_face_labeling(oldmodel),tags="interior")

#using Gridap.Visualization
#writevtk(oldtrian,"oldtrian")
#writevtk(btrian,"btrian",cellfields=["normal"=>nb],celldata=["oldcell"=>get_cell_id(btrian)])
#writevtk(strian,"strian",cellfields=["normal"=>ns],
#  celldata=["oldcell_left"=>get_cell_id(strian).plus,"oldcell_right"=>get_cell_id(strian).minus])

end # module
