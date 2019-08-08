module CartesianGeometryTests

using Test
using Gridap
using Gridap.DiscreteModels: dict_to_model
using JSON

grid = CartesianGrid(domain=(0.0,1.0,-1.0,2.0),partition=(3,4))
test_grid(grid,20,12)

@test celldim(grid) == 2
@test pointdim(grid) == 2
@test ncells(grid) == 3*4
@test npoints(grid) == (3+1)*(4+1)

x = points(grid)

@test isa(x,CellValue{Point{2,Float64}})

@test length(x) == 20
@test x[13] == Point(0.0,1.25)
@test x[4] == Point(1.0,-1.0)
@test x[1] == Point(0.0,-1.0)
@test x[end] == Point(1.0,2.0)

t = cells(grid)

@test isa(t,CellArray{Int,1})

@test length(t) == 12

@test t[1] == [1,2,5,6]
@test t[11] == [14,15,18,19]

c = celltypes(grid)
@test isa(c,ConstantCellValue{NTuple{2,Int}})
@test c.value == (HEX_AXIS,HEX_AXIS)
@test length(c) == 12

o = cellorders(grid)
@test isa(o,ConstantCellValue{Int})
@test o.value == 1
@test length(o) == 12

trian = Triangulation(grid)

xe = CellPoints(trian)
@test isa(xe,CellPoints{2,Float64})

cb = CellBasis(trian)
@test isa(cb,CellBasis{2,Float64})

model = CartesianDiscreteModel(
  domain=(0.0,1.0,-1.0,2.0,0.0,1.0),
  partition=(3,4,2))

test_discrete_model(model,3)

grid3 = Grid(model,3)
np = (3+1)*(4+1)*(2+1)
nc = 3*4*2
test_grid(grid3,np,nc)

@test isa(grid3,CartesianGrid{3})

gridgraph = FullGridGraph(model)

@test isa(gridgraph, FullGridGraph)

grid2 = Grid(model,2)
grid1 = Grid(model,1)
grid0 = Grid(model,0)

labels = FaceLabels(model)

@test isa(labels,FaceLabels)
@test length(labels_on_dim(labels,3)) == ncells(grid3)
@test length(labels_on_dim(labels,2)) == ncells(grid2)
@test length(labels_on_dim(labels,1)) == ncells(grid1)
@test length(labels_on_dim(labels,0)) == ncells(grid0)

@test ntags(labels) == 28
@test name_from_tag(labels,28) == "boundary"

cgrid = CartesianGrid(domain=(0.0,1.0,-1.0,2.0),partition=(3,4))

grid = FlexibleUnstructuredGrid(cgrid)
test_grid(grid,20,12)

grid = UnstructuredGrid(cgrid)
test_grid(grid,20,12)

# Serialization

model = CartesianDiscreteModel(partition=(2,2,2))
s = json(model)

dict = JSON.parse(s)
model = dict_to_model(dict)
test_discrete_model(model,3)

d = mktempdir()
filename = joinpath(d,"model.json")

open(filename,"w") do f
  JSON.print(f,model)
end

model = DiscreteModelFromFile(filename)
test_discrete_model(model,3)

end # module
