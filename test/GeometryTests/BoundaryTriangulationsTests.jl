module BoundaryTriangulationsTests

using Test
using Gridap.Helpers
using Gridap.Fields
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using FillArrays

domain = (0,4,0,4)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

btrian = BoundaryTriangulation(model)

face_s_q = get_cell_ref_map(btrian)

s1 = Point(0.0)
s2 = Point(0.5)
s = [s1,s2]
face_to_s = Fill(s,length(face_s_q))

face_to_q = lazy_map(evaluate,face_s_q,face_to_s)
@test isa(face_to_q,Geometry.FaceCompressedVector)

cell_shapefuns = get_cell_shapefuns(btrian.cell_trian)

face_shapefuns = lazy_map(Reindex(cell_shapefuns),get_cell_id(btrian))

face_shapefuns_q = lazy_map(evaluate,face_shapefuns,face_to_q)
@test isa(face_shapefuns_q,Geometry.FaceCompressedVector)

face_to_nvec = get_facet_normal(btrian)

face_to_nvec_s = lazy_map(evaluate,face_to_nvec,face_to_s)

print_op_tree(face_to_nvec_s)

display(face_to_nvec_s)


#print_op_tree(face_shapefuns)

#print_op_tree(face_to_q)
#print_op_tree(face_shapefuns_q)

kk


test_triangulation(btrian)

kk

cellids = collect(1:num_cells(model))

face_to_cellid = reindex(cellids,btrian)
@test face_to_cellid == get_face_to_cell(btrian)

trian = get_volume_triangulation(btrian)
q2x = get_cell_map(trian)
s2x = restrict(q2x,btrian)

s = CompressedArray([Point{1,Float64}[(0.25,),(0.75,)]],get_cell_type(btrian))
x = evaluate(s2x,s)
r = Vector{Point{2,Float64}}[
  [(0.5,0.0),(1.5,0.0)],[(0.0,0.5),(0.0,1.5)],
  [(2.5,0.0),(3.5,0.0)],[(4.0,0.5),(4.0,1.5)],
  [(0.5,4.0),(1.5,4.0)],[(0.0,2.5),(0.0,3.5)],
  [(2.5,4.0),(3.5,4.0)],[(4.0,2.5),(4.0,3.5)]]
test_array(x,r)

nvec = get_normal_vector(btrian)
nvec_x = evaluate(nvec,s)

domain = (0,4,0,4)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

btrian = BoundaryTriangulation(model,"tag_8")
test_boundary_triangulation(btrian)
@test btrian.face_trian.cell_to_oldcell == get_face_to_face(btrian)
@test 1 == get_cell_around(btrian)

btrian = BoundaryTriangulation(model,get_face_labeling(model),"tag_8")
test_boundary_triangulation(btrian)

s = CompressedArray([Point{1,Float64}[(0.25,),(0.75,)]],get_cell_type(btrian))
nvec = get_normal_vector(btrian)
nvec_x = evaluate(nvec,s)
s2x = get_cell_map(btrian)
x = evaluate(s2x,s)

@test get_cell_id(btrian) == get_face_to_cell(btrian)
r = rand(num_cells(trian))
@test reindex(r,btrian) == r[get_face_to_cell(btrian)]

oldbtrian = BoundaryTriangulation(oldmodel)

bface_to_oldbface = collect(1:32)
btrian = TriangulationPortion(oldbtrian,bface_to_oldbface)
test_triangulation(btrian)

nb = get_normal_vector(btrian)
@test length(nb) == num_cells(btrian)


end # module
