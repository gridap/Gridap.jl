module WrappersTests

using Test
using Numa
using Numa.FieldValues
using Numa.CellValues

include("Helpers.jl")

v = [ Point{2}(i,i) for i in 1:10 ]
cv = CellValueFromArray(v)
test_index_cell_value( cv, v )

v = CartesianIndices((3,10))
cv = CellValueFromArray(v)
test_index_cell_value( cv, v )

v = [ [i,i,i] for i in 1:10 ]
cv = CellValueFromArray(v)
test_index_cell_array( cv, v )



end # module WrappersTests


#@testset "Wrappers" begin
#
#
#  a = [ [i+1,i+2,i+3] for i in 1:10 ]
#
#  ca = CellArrayFromArrayOfArrays(a)
#
#  @test length(ca) == length(a)
#  @test size(ca) == size(a)
#  for (cai,ai) in zip(ca,a)
#    @assert cai === ai
#  end
#  @test IndexStyle(typeof(a)) == IndexStyle(typeof(ca))
#
#  #       1,2,3,4,5,6,7,8,9,0,1,2
#  data = [2,3,1,3,6,7,3,2,5,6,3,4]
#  ptrs = [1,4,4,7,13]
#
#  ca = CellVectorFromDataAndPtrs(data,ptrs)
#
#  @test length(ca) == length(ptrs)-1
#  @test ca[1] == data[1:3]
#  @test ca[2] == data[4:3]
#  @test ca[3] == data[4:6]
#  @test ca[4] == data[7:12]
#
#  @test cellsize(ca) == (6,)
#
#  ca = CellVectorFromDataAndStride(data,3)
#
#  @test length(ca) == 4
#  @test ca[1] == data[1:3]
#  @test ca[2] == data[4:6]
#  @test ca[3] == data[7:9]
#  @test ca[4] == data[10:12]
#
#  @test cellsize(ca) == (3,)
#
#  p = Point(2.0,1.0)
#  #       1,2,3,4,5,6,7,8,9,0,1,2
#  data = [2,3,1,3,4,4,3,2,5,4,3,4]
#  ptrs = [1,4,4,7,13]
#  gid_to_val = [p, 2*p, 3*p, -p, p]
#
#  lid_to_gid = CellVectorFromDataAndPtrs(data,ptrs)
#  ca = CellVectorFromLocalToGlobal(lid_to_gid,gid_to_val)
#
#  @test length(ca) == length(ptrs)-1
#  @test cellsize(ca) == (6,)
#
#  @test ca[1] == gid_to_val[data[1:3]]
#  @test ca[2] == gid_to_val[data[4:3]]
#  @test ca[3] == gid_to_val[data[4:6]]
#  @test ca[4] == gid_to_val[data[7:12]]
#
#  #       1,2,3,4,5,6,7,8,9,0,1,2
#  data = [-2,3,-1,3,4,-4,-3,2,5,4,3,-4]
#  ptrs = [1,4,4,7,13]
#  p = 1
#  gid_to_val_pos = [p, 2*p, 3*p, -p, p]
#  gid_to_val_neg = [2*p, p, 5*p, 3*p, p]
#
#  lid_to_gid = CellVectorFromDataAndPtrs(data,ptrs)
#  ca = CellVectorFromLocalToGlobalPosAndNeg(lid_to_gid,gid_to_val_pos,gid_to_val_neg)
#
#  @test length(ca) == length(ptrs)-1
#  @test cellsize(ca) == (6,)
#
#  @test ca[1] == [1, 3, 2]
#  @test ca[2] == Int64[]
#  @test ca[3] == [3, -1, 3]
#  @test ca[4] == [5, 2, 1, -1, 3, 3]
#
#  cell_to_x_l = [2,3,1,3,4,4,3,2,5,4,3,4]
#  cell_to_x_p = [1,4,4,7,13]
#  cell_to_x = CellVectorFromDataAndPtrs(cell_to_x_l,cell_to_x_p)
#  x_to_vals_l = [5,4,1,2,3,6,7,8,9,10]
#  x_to_vals_p = [1,3,5,6,10,11]
#  x_to_vals = CellVectorFromDataAndPtrs(x_to_vals_l,x_to_vals_p)
#  using Numa.CellValues: CellVectorByComposition
#  cell_to_vals = CellVectorByComposition(cell_to_x, x_to_vals)
#  @test cell_to_vals[1] == [1,2,3,5,4]
#  
#end
