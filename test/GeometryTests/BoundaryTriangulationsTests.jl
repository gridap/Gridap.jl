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
model = simplexify(CartesianDiscreteModel(domain,partition))

btrian = BoundaryTriangulation(model,tags=[7,8])
test_triangulation(btrian)
@test get_background_model(btrian) === model
glue = get_glue(btrian,Val(1))
@test glue.tface_to_mface === btrian.glue.face_to_bgface
glue = get_glue(btrian,Val(2))
@test glue.tface_to_mface === btrian.glue.face_to_cell

face_s_q = glue.tface_to_mface_map

s1 = Point(0.0)
s2 = Point(0.5)
s = [s1,s2]
face_to_s = Fill(s,length(face_s_q))

face_to_q = lazy_map(evaluate,face_s_q,face_to_s)
@test isa(face_to_q,Geometry.FaceCompressedVector)

cell_grid = get_grid(get_background_model(btrian))
cell_shapefuns = get_cell_shapefuns(cell_grid)
cell_grad_shapefuns = lazy_map(Broadcasting(∇),cell_shapefuns)

face_shapefuns = lazy_map(Reindex(cell_shapefuns),glue.tface_to_mface)
face_grad_shapefuns = lazy_map(Reindex(cell_grad_shapefuns),glue.tface_to_mface)

face_shapefuns_q = lazy_map(evaluate,face_shapefuns,face_to_q)
test_array(face_shapefuns_q,collect(face_shapefuns_q))
@test isa(face_shapefuns_q,Geometry.FaceCompressedVector)

face_grad_shapefuns_q = lazy_map(evaluate,face_grad_shapefuns,face_to_q)
test_array(face_grad_shapefuns_q,collect(face_grad_shapefuns_q))
@test isa(face_grad_shapefuns_q,Geometry.FaceCompressedVector)

face_to_nvec = get_facet_normal(btrian)
face_to_nvec_s = lazy_map(evaluate,face_to_nvec,face_to_s)
test_array(face_to_nvec_s,collect(face_to_nvec_s))

##print_op_tree(face_shapefuns)
#print_op_tree(face_to_q)
#print_op_tree(face_shapefuns_q)
#print_op_tree(face_grad_shapefuns_q)
#print_op_tree(face_to_nvec_s)

Ω = Triangulation(model)
Γ = btrian
∂Γ = BoundaryTriangulation(Γ)
test_triangulation(∂Γ)
@test ∂Γ.rtrian === Γ
@test isa(∂Γ.dtrian,BoundaryTriangulation)
glue = get_glue(∂Γ,Val(2))
glue = get_glue(∂Γ,Val(1))
@test is_change_possible(Ω,∂Γ)
@test is_change_possible(Γ,∂Γ)

domain = (0,4,0,4)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

btrian = BoundaryTriangulation(model)
test_triangulation(btrian)
@test get_background_model(btrian) === model
glue = get_glue(btrian,Val(1))
@test glue.tface_to_mface === btrian.glue.face_to_bgface
glue = get_glue(btrian,Val(2))
@test glue.tface_to_mface === btrian.glue.face_to_cell

@test isa(get_cell_map(btrian)[1],AffineField)

face_s_q = glue.tface_to_mface_map

s1 = Point(0.0)
s2 = Point(0.5)
s = [s1,s2]
face_to_s = Fill(s,length(face_s_q))

face_to_q = lazy_map(evaluate,face_s_q,face_to_s)
@test isa(face_to_q,Geometry.FaceCompressedVector)

cell_grid = get_grid(get_background_model(btrian))
cell_shapefuns = get_cell_shapefuns(cell_grid)
cell_grad_shapefuns = lazy_map(Broadcasting(∇),cell_shapefuns)

face_shapefuns = lazy_map(Reindex(cell_shapefuns),glue.tface_to_mface)
face_grad_shapefuns = lazy_map(Reindex(cell_grad_shapefuns),glue.tface_to_mface)

face_shapefuns_q = lazy_map(evaluate,face_shapefuns,face_to_q)
test_array(face_shapefuns_q,collect(face_shapefuns_q))
@test isa(face_shapefuns_q,Geometry.FaceCompressedVector)

face_grad_shapefuns_q = lazy_map(evaluate,face_grad_shapefuns,face_to_q)
test_array(face_grad_shapefuns_q,collect(face_grad_shapefuns_q))
@test isa(face_grad_shapefuns_q,Geometry.FaceCompressedVector)

face_to_nvec = get_facet_normal(btrian)
face_to_nvec_s = lazy_map(evaluate,face_to_nvec,face_to_s)
test_array(face_to_nvec_s,collect(face_to_nvec_s))
@test isa(face_to_nvec_s,Geometry.FaceCompressedVector)

face_to_tvec = get_edge_tangent(btrian)
face_to_tvec_s = lazy_map(evaluate,face_to_tvec,face_to_s)
test_array(face_to_tvec_s,collect(face_to_tvec_s))
@test isa(face_to_tvec_s,Geometry.FaceCompressedVector)

#print_op_tree(face_shapefuns_q)
#print_op_tree(face_grad_shapefuns_q)
#print_op_tree(face_to_nvec_s)

s = CompressedArray([Point{1,Float64}[(0.25,),(0.75,)]],get_cell_type(btrian))

s2q = glue.tface_to_mface_map

q = lazy_map(evaluate,s2q,s)
r = Vector{Point{2,Float64}}[
  [(0.25,0.0),(0.75,0.0)],[(0.0,0.25),(0.0,0.75)],
  [(0.25,0.0),(0.75,0.0)],[(1.0,0.25),(1.0,0.75)],
  [(0.25,1.0),(0.75,1.0)],[(0.0,0.25),(0.0,0.75)],
  [(0.25,1.0),(0.75,1.0)],[(1.0,0.25),(1.0,0.75)]]
test_array(q,r)

q2x = get_cell_map(cell_grid)
s2x = lazy_map(∘,lazy_map(Reindex(q2x),glue.tface_to_mface),s2q)

x = lazy_map(evaluate,s2x,s)
r = Vector{Point{2,Float64}}[
  [(0.5,0.0),(1.5,0.0)],[(0.0,0.5),(0.0,1.5)],
  [(2.5,0.0),(3.5,0.0)],[(4.0,0.5),(4.0,1.5)],
  [(0.5,4.0),(1.5,4.0)],[(0.0,2.5),(0.0,3.5)],
  [(2.5,4.0),(3.5,4.0)],[(4.0,2.5),(4.0,3.5)]]
test_array(x,r)

nvec = get_facet_normal(btrian)

nvec_s = lazy_map(evaluate,nvec,s)
r = Vector{Point{2,Float64}}[
  [(0.0,-1.0),(0.0,-1.0)],[(-1.0,0.0),(-1.0,0.0)],
  [(0.0,-1.0),(0.0,-1.0)],[(1.0,-0.0),(1.0,-0.0)],
  [(0.0,1.0),(0.0,1.0)],[(-1.0,0.0),(-1.0,0.0)],
  [(0.0,1.0),(0.0,1.0)],[(1.0,-0.0),(1.0,-0.0)]]
test_array(nvec_s,r)

tvec = get_edge_tangent(btrian)

tvec_s = lazy_map(evaluate,tvec,s)

rotate_90(x) = [Point(x[1][2], -x[1][1]), Point(x[2][2], -x[2][1])]

rr = rotate_90.(r)

test_array(rr, tvec_s)

cellids = collect(1:num_cells(model))

face_to_cellid = lazy_map(Reindex(cellids),glue.tface_to_mface)
@test face_to_cellid == glue.tface_to_mface

q2x = get_cell_map(cell_grid)
s2q = glue.tface_to_mface_map
s2x = lazy_map(∘,lazy_map(Reindex(q2x),glue.tface_to_mface),s2q)

s = CompressedArray([Point{1,Float64}[(0.25,),(0.75,)]],get_cell_type(btrian))
x = lazy_map(evaluate,s2x,s)
r = Vector{Point{2,Float64}}[
  [(0.5,0.0),(1.5,0.0)],[(0.0,0.5),(0.0,1.5)],
  [(2.5,0.0),(3.5,0.0)],[(4.0,0.5),(4.0,1.5)],
  [(0.5,4.0),(1.5,4.0)],[(0.0,2.5),(0.0,3.5)],
  [(2.5,4.0),(3.5,4.0)],[(4.0,2.5),(4.0,3.5)]]
test_array(x,r)

nvec = get_facet_normal(btrian)
nvec_x = lazy_map(evaluate,nvec,s)

domain = (0,4,0,4)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

btrian = BoundaryTriangulation(model,tags="tag_8")
test_triangulation(btrian)
@test all(1 .== btrian.glue.face_to_lcell)

btrian = BoundaryTriangulation(model,get_face_labeling(model),tags="tag_8")
test_triangulation(btrian)

s = CompressedArray([Point{1,Float64}[(0.25,),(0.75,)]],get_cell_type(btrian))
nvec = get_facet_normal(btrian)
nvec_x = lazy_map(evaluate,nvec,s)
s2x = get_cell_map(btrian)
x = lazy_map(evaluate,s2x,s)

btrian = BoundaryTriangulation(model,tags=["tag_6","tag_7"])
test_triangulation(btrian)

ids = [1,3,4]
vtrian = view(btrian,ids)
test_triangulation(vtrian)
@test num_cells(vtrian) == 3
bglue = get_glue(btrian,Val(1))
vglue = get_glue(vtrian,Val(1))
@test vglue.tface_to_mface == bglue.tface_to_mface[ids]
bglue = get_glue(btrian,Val(2))
vglue = get_glue(vtrian,Val(2))
@test vglue.tface_to_mface == bglue.tface_to_mface[ids]

end # module
