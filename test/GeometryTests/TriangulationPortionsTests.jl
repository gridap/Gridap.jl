module TriangulationPortionsTests

using Gridap.ReferenceFEs
using Gridap.Geometry

domain = (0,1,0,1)
partition = (10,10)
oldgrid = CartesianGrid(domain,partition)

cell_to_oldcell = collect(1:34)
trian = TriangulationPortion(oldgrid,cell_to_oldcell)
test_triangulation(trian)

end # module

