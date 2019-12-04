module VtkTests

using Gridap.Geometry: ConformingTrianMock
using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Visualization

d = mktempdir()
f = joinpath(d,"trian")

trian = ConformingTrianMock()

node_ids = collect(1:num_nodes(trian))
cell_ids = collect(1:num_cells(trian))

mean(x) = sum(x)/length(x)

cell_center = apply(mean, get_cell_coordinates(trian) )

write_vtk_file(trian,f,nodaldata=["nodeid"=>node_ids],celldata=["cellid"=>cell_ids,"centers"=>cell_center])

reffe = LagrangianRefFE(VectorValue{3,Float64},WEDGE,(3,3,4))
f = joinpath(d,"reffe")
writevtk(reffe,f)

f = joinpath(d,"poly")
writevtk(HEX,f)

rm(d,recursive=true)

end # module
