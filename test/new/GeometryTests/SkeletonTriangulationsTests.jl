module SkeletonTriangulationsTests

using Test
using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry

#function polar(q)
#  r, t, z = q
#  x = r*cos(t)
#  y = r*sin(t)
#  Point(x,y,z)
#end
#domain = (1,2,0,pi,0,0.5)
#partition = (10,30,4)
#model = CartesianDiscreteModel(domain,partition,polar)

domain = (0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)

strian = SkeletonTriangulation(model)
test_triangulation(strian)

trian = get_volume_triangulation(strian)

s = CompressedArray([Point{2,Float64}[(0.25,0.25),(0.75,0.75)]],get_cell_type(strian))

s2x = get_cell_map(strian)
x = evaluate(s2x,s)

nvec = get_normal_vector(strian)
nx = evaluate(nvec,s)
collect(nx)

fun(x) = sin(pi*x[1])*cos(pi*x[2])

q2x = get_cell_map(trian)

funq = compose(fun,q2x)

fun_gamma = restrict(funq,strian)

@test isa(fun_gamma, SkeletonPair)

cellids = collect(1:num_cells(trian))

cellids_gamma = reindex(cellids,strian)
@test isa(fun_gamma, SkeletonPair)
@test cellids_gamma.left == get_face_to_cell(strian.left)
@test cellids_gamma.right == get_face_to_cell(strian.right)

#using Gridap.Visualization
#
#writevtk(trian,"trian")
#writevtk(strian,"strian")
#writevtk(x,"x",nodaldata=["nvec" => nx])

end # module
