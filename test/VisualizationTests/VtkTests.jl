module VtkTests

using Test
using Gridap.Geometry: GridMock
using Gridap.Geometry: DiscreteModelMock
using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Visualization

d = mktempdir()
f = joinpath(d,"trian")

trian = GridMock()

node_ids = collect(1:num_nodes(trian))
cell_ids = collect(1:num_cells(trian))

mean(x) = sum(x)/length(x)

cell_center = apply(mean, get_cell_coordinates(trian) )

write_vtk_file(
  trian,f,
  nodaldata=["nodeid"=>node_ids],
  celldata=["cellid"=>cell_ids,"centers"=>cell_center])

reffe = LagrangianRefFE(VectorValue{3,Float64},WEDGE,(3,3,4))
f = joinpath(d,"reffe")
writevtk(reffe,f)

f = joinpath(d,"poly")
writevtk(HEX,f)

domain = (0,1,0,1,0,1)
partition = (3,4,2)
grid = CartesianGrid(domain,partition)
model = UnstructuredDiscreteModel(grid)

f = joinpath(d,"model")
writevtk(model,f)

domain = (0,1,0,1,0,1)
partition = (3,4,2)
model = CartesianDiscreteModel(domain,partition)

f = joinpath(d,"model")
writevtk(model,f)

f = joinpath(d,"model")
model = DiscreteModelMock()
writevtk(model,get_face_labeling(model),f)

f = joinpath(d,"trian")
trian = GridMock()
writevtk(trian,f,order=2)

domain = (0,1,0,1)
partition = (2,4)
trian = CartesianGrid(domain,partition)

writevtk(trian,f,nsubcells=5,celldata=["rnd"=>rand(num_cells(trian))])

fun(x) = sin(4*x[1]*pi)*cos(5*x[2]*pi)
cf = compose(fun, get_cell_map(trian))

writevtk(trian,f,nsubcells=10, cellfields=["cf" => cf])

trian = GridMock()

p1 = Point{2,Float64}[(0.25,0.25),(0.75,0.75)]
p2 = Point{2,Float64}[(0.2,0.2),(0.4,0.4)]
q = CompressedArray([p1,p2], get_cell_type(trian))
q2x = get_cell_map(trian)
x = evaluate(q2x,q)

f = joinpath(d,"x")
writevtk(x,f,celldata=["cellid" => collect(1:num_cells(trian))], nodaldata = ["x" => x])

f = joinpath(d,"trian")
writevtk(trian,f)

rm(d,recursive=true)

end # module
