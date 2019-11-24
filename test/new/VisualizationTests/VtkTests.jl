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

#reffe = LagrangianRefFE(Float64,QUAD,2)
#
#"""
#"""
#function writevtk(reffe::NodalReferenceFE, filebase)
#
#  grid = UnstructuredGrid(reffe)
#
#  face_to_own_nodes = get_face_own_nodeids(reffe)
#  node_to_owner = zeros(Int,num_nodes(reffe))
#  for (face, nodeids) in enumerate(face_to_own_nodes)
#    node_to_owner[nodeids] .= face
#  end
#
#  writevtk_file(grid,filebase;nodaldata=["face_owner"=>node_to_owner])
#
#end
#
#writevtk(reffe,"reffe")

end # module
