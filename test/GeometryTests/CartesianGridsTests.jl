module CartesianGridsTests

using FillArrays
using Test

using Gridap.Helpers
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.TensorValues: det
using Gridap.Fields
using Gridap.Geometry
using Gridap.ReferenceFEs

domain = (0,1)
partition = (1,)
grid = CartesianGrid(domain,partition)
test_grid(grid)
test_grid(UnstructuredGrid(grid))

pmin = Point(0.0,-1.0)
pmax = Point(1.0,2.0)
partition = (3,4)
grid1 = CartesianGrid(pmin,pmax,partition)

domain = (0.0,1.0,-1.0,2.0)
partition = (3,4)
grid2 = CartesianGrid(domain,partition)

domain = [0.0,1.0,-1.0,2.0]
partition = [3,4]
grid3 = CartesianGrid(domain,partition)

origin = Point(0,-1)
sizes = (1/3,3/4)
partition = (3,4)
grid4 = CartesianGrid(origin,sizes,partition)

for grid in (grid1,grid2,grid3,grid4)
  test_grid(grid)
  @test is_oriented(grid) == true

  desc = get_cartesian_descriptor(grid)
  @test desc.origin == origin
  @test desc.sizes == sizes
  @test desc.partition == partition

  @test num_cell_dims(grid) == 2
  @test num_point_dims(grid) == 2
  @test num_cells(grid) == 3*4
  @test num_nodes(grid) == (3+1)*(4+1)

  x = get_node_coordinates(grid)
  @test length(x) == 20
  @test x[13] == Point(0.0,1.25)
  @test x[4] == Point(1.0,-1.0)
  @test x[1] == Point(0.0,-1.0)
  @test x[end] == Point(1.0,2.0)

  t = get_cell_node_ids(grid)
  @test t[1] == [1,2,5,6]
  @test t[11] == [14,15,18,19]

  cell_type = get_cell_type(grid)
  @test isa(cell_type,Fill)

  reffes = get_reffes(grid)
  @test length(reffes) == 1
  @test reffes[1] == QUAD4

  ugrid = UnstructuredGrid(grid)
  test_grid(ugrid)
  @test is_oriented(ugrid) == true
end

domain = (0,1,0,1)
partition = (2,2)

grid = CartesianGrid(domain,partition)

map = get_cell_map(grid)
∇map = lazy_map(∇,map)

x = [Point(0.5,0.5),]
ax = Fill(x,partition)
r = Vector{Point{2,Float64}}[[(0.25, 0.25)], [(0.75, 0.25)], [(0.25, 0.75)], [(0.75, 0.75)]]
∇r = Vector{TensorValue{2,2,Float64,4}}[
  [(0.5, 0.0, 0.0, 0.5)], [(0.5, 0.0, 0.0, 0.5)],
  [(0.5, 0.0, 0.0, 0.5)], [(0.5, 0.0, 0.0, 0.5)]]
#test_array_of_fields(map,ax,r,grad=∇r)
@test all(lazy_map(test_field,map,ax,r))
@test all(lazy_map(test_field,∇map,ax,∇r))

mx = lazy_map(evaluate,map,ax)
∇mx = lazy_map(evaluate,∇map,ax)

test_array(mx,reshape(r,size(map)))
test_array(∇mx,reshape(∇r,size(map)))
@test isa(∇mx,Fill)

ptgrid = simplexify(grid,positive=true)
test_grid(ptgrid)

Jtk = lazy_map(∇,get_cell_map(ptgrid))
X0k = lazy_map(first,get_cell_coordinates(ptgrid))
Jt0k = lazy_map(evaluate,Jtk,X0k)
@test all(>(0),lazy_map(det,Jt0k))

domain = (0,1,0,1,0,1)
partition = (2,2,2)
grid = CartesianGrid(domain,partition)
tgrid = simplexify(grid)
test_grid(tgrid)

ptgrid = simplexify(grid,positive=true)
test_grid(ptgrid)

Jtk = lazy_map(∇,get_cell_map(ptgrid))
X0k = lazy_map(first,get_cell_coordinates(ptgrid))
Jt0k = lazy_map(evaluate,Jtk,X0k)
@test all(>(0),lazy_map(det,Jt0k))

grid = compute_reference_grid(HEX8,4)
test_grid(grid)

grid = compute_linear_grid(HEX8)
test_grid(grid)

end # module
