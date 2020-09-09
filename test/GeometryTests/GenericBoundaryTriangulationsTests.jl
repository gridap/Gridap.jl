module GenericBoundaryTriangulationsTests

using Test
using Gridap.TensorValues
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using FillArrays

using Gridap.Geometry: FaceCellCoordinates
using Gridap.Geometry: FaceToCellMapValue
using Gridap.Geometry: CellShapeFunsAtFaces

domain = (0,4,0,4)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
model = UnstructuredDiscreteModel(model)
topo = get_grid_topology(model)

btrian = GenericBoundaryTriangulation(model,get_isboundary_face(topo,1))
test_triangulation(btrian)

s = CompressedArray([Point{1,Float64}[(0.25,),(0.75,)]],get_cell_type(btrian))

face_to_fvertex_to_qcoors = FaceCellCoordinates(btrian)
r =  Vector{Point{2,Float64}}[
  [(0.0,0.0),(1.0,0.0)],[(0.0,0.0),(0.0,1.0)],
  [(0.0,0.0),(1.0,0.0)],[(1.0,0.0),(1.0,1.0)],
  [(0.0,1.0),(1.0,1.0)],[(0.0,0.0),(0.0,1.0)],
  [(0.0,1.0),(1.0,1.0)],[(1.0,0.0),(1.0,1.0)]]
test_array(face_to_fvertex_to_qcoors,r)

s2q = get_face_to_cell_map(btrian)

q = evaluate(s2q,s)
r = Vector{Point{2,Float64}}[
  [(0.25,0.0),(0.75,0.0)],[(0.0,0.25),(0.0,0.75)],
  [(0.25,0.0),(0.75,0.0)],[(1.0,0.25),(1.0,0.75)],
  [(0.25,1.0),(0.75,1.0)],[(0.0,0.25),(0.0,0.75)],
  [(0.25,1.0),(0.75,1.0)],[(1.0,0.25),(1.0,0.75)]]
test_array(q,r)
@test isa(q,FaceToCellMapValue)

#q2x = reindex(get_cell_map(get_grid(model)), btrian.glue.face_to_cell)
q2x = get_cell_map(get_grid(model))

s2x = q2x∘s2q

x = evaluate(s2x,s)
r = Vector{Point{2,Float64}}[
  [(0.5,0.0),(1.5,0.0)],[(0.0,0.5),(0.0,1.5)],
  [(2.5,0.0),(3.5,0.0)],[(4.0,0.5),(4.0,1.5)],
  [(0.5,4.0),(1.5,4.0)],[(0.0,2.5),(0.0,3.5)],
  [(2.5,4.0),(3.5,4.0)],[(4.0,2.5),(4.0,3.5)]]
test_array(x,r)

cell_shapefuns_q = reindex(get_cell_shapefuns(get_grid(model)), btrian.glue.face_to_cell)
cell_shapefuns_s = cell_shapefuns_q∘s2q
@test isa(evaluate(cell_shapefuns_s,s), CellShapeFunsAtFaces)

nvec = get_normal_vector(btrian)

nvec_s = evaluate(nvec∘get_cell_map(btrian),s)
r = Vector{Point{2,Float64}}[
  [(0.0,-1.0),(0.0,-1.0)],[(-1.0,0.0),(-1.0,0.0)],
  [(0.0,-1.0),(0.0,-1.0)],[(1.0,-0.0),(1.0,-0.0)],
  [(0.0,1.0),(0.0,1.0)],[(-1.0,0.0),(-1.0,0.0)],
  [(0.0,1.0),(0.0,1.0)],[(1.0,-0.0),(1.0,-0.0)]]
test_array(nvec_s,r)

end # module
