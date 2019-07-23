module GridGraphsTests

using Test
using Gridap
using Gridap.CellValuesGallery
using Gridap.GridGraphs: GridGraphFromData
using Gridap.GridGraphs: FullGridGraphFromData

c2v = [[1,2,4,5],[2,3,5,6]]
v2c = [[1,],[1,2],[2,],[1,],[1,2],[2,]]
c2e = [[1,2,3,4],[5,6,4,7]]
e2c = [[1,],[1,],[1,],[1,2],[2,],[2,],[2,]]
e2v = [[1,2],[4,5],[1,4],[2,5],[2,3],[5,6],[3,7]]
v2e = [[1,3],[1,4,5],[5,7],[2,3],[2,4,6],[6,7]]
v2v = [[1,],[2,],[3,],[4,],[5,],[6,]]
e2e = [[1,],[2,],[3,],[4,],[5,],[6,],[7,]]
c2c = [[1,],[2,]]

c2v = CellValueFromArray(c2v)
v2c = CellValueFromArray(v2c)
c2e = CellValueFromArray(c2e)
e2c = CellValueFromArray(e2c)
e2v = CellValueFromArray(e2v)
v2e = CellValueFromArray(v2e)
v2v = CellValueFromArray(v2v)
e2e = CellValueFromArray(e2e)
c2c = CellValueFromArray(c2c)

primal = [c2v,c2e,c2c]
dual = [v2c,e2c,c2c]
graph = GridGraphFromData(primal,dual)
nd = 2
test_grid_graph(graph,nd)

data = Matrix{IndexCellVector}(undef,(3,3))
data[1,:] = [v2v,v2e,v2c]
data[2,:] = [e2v,e2e,e2c]
data[3,:] = [c2v,c2e,c2c]
graph = FullGridGraphFromData(data)
test_full_grid_graph(graph,nd)

end # module
