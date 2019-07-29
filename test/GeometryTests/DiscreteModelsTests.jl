module DiscreteModelsTests

using Test
using Gridap
using Gridap.DiscreteModels: DiscreteModelFromData
using Gridap.CellValuesGallery
using Gridap.Grids: GridFromData
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
v2x = Point{2,Float64}[(0.,0.),(0.,1.),(0.,2.),(1.,0.),(1.,1.),(1.,2.)]

c2v = CellValueFromArray(c2v)
v2c = CellValueFromArray(v2c)
c2e = CellValueFromArray(c2e)
e2c = CellValueFromArray(e2c)
e2v = CellValueFromArray(e2v)
v2e = CellValueFromArray(v2e)
v2v = CellValueFromArray(v2v)
e2e = CellValueFromArray(e2e)
c2c = CellValueFromArray(c2c)
v2x = CellValueFromArray(v2x)

data = Matrix{IndexCellVector}(undef,(3,3))
data[1,:] = [v2v,v2e,v2c]
data[2,:] = [e2v,e2e,e2c]
data[3,:] = [c2v,c2e,c2c]
graph = FullGridGraphFromData(data)

order = 1

nc = length(c2v)
t2 = (HEX_AXIS,HEX_AXIS)
c2t = ConstantCellValue(t2,nc)
c2o = ConstantCellValue(order,nc)
grid2 = GridFromData(v2x,c2v,c2t,c2o)

ne = length(e2v)
t1 = (HEX_AXIS,)
e2t = ConstantCellValue(t1,ne)
e2o = ConstantCellValue(order,ne)
grid1 = GridFromData(v2x,e2v,e2t,e2o)

nv = length(v2v)
t0 = ()
v2t = ConstantCellValue(t0,nv)
v2o = ConstantCellValue(order,nv)
grid0 = GridFromData(v2x,v2v,v2t,v2o)

v2l = CellValueFromArray([1,1,1,2,2,1])
e2l = CellValueFromArray([2,2,3,3,2,2,1])
c2l = CellValueFromArray([5,5])

dim_to_labels = [v2l, e2l, c2l]
tag1 = [1,3]
tag2 = [5,3,2]
tags = [tag1,tag2]
tag_to_name = ["tag1","tag2"]

facelabels = FaceLabels(dim_to_labels,tags,tag_to_name)

grids = Grid[grid0,grid1,grid2]

model = DiscreteModelFromData(grids,graph,facelabels)

test_discrete_model(model,2)

vgrid = Grid(model)

end # module
