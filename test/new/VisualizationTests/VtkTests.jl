module VtkTests

using Gridap.Geometry: ConformingTrianMock
using Gridap.Arrays
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

writevtk_file(trian,f,nodaldata=["nodeid"=>node_ids],celldata=["cellid"=>cell_ids,"centers"=>cell_center])

rm(d,recursive=true)

end # module
