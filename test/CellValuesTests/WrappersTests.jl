module WrappersTests

using Test
using Gridap
using Gridap.FieldValues
using Gridap.CellValues
using Gridap.CellValues.Testers
using Gridap.CellValues.Wrappers

include("CartesianArray.jl")

v = [ Point{2}(i,i) for i in 1:10 ]
cv = CellValueFromArray(v)
test_index_cell_value( cv, v )

v = CartesianIndices((3,10))
cv = CellValueFromArray(v)
test_index_cell_value( cv, v )

v = [ [i,i,i] for i in 1:10 ]
cv = CellValueFromArray(v)
test_index_cell_array( cv, v )

r = [[2, 2, 2], [4, 4, 4], [6, 6, 6], [8, 8, 8],
     [10, 10, 10], [12, 12, 12], [14, 14, 14],
     [16, 16, 16], [18, 18, 18], [20, 20, 20]]
@test r == collect(cv+cv)

a = [ [i+1,i+2,i+3] for i in 1:10 ]
ca = CellValueFromArray(a)
test_index_cell_array( ca, a )

data = [2,3,1,3,6,7,3,2,5,6,3,4]
ptrs = [1,4,4,7,13]
ca = CellVectorFromDataAndPtrs(data,ptrs)
a = [ data[ptrs[i]:ptrs[i+1]-1] for i in 1:length(ptrs)-1]
test_index_cell_array( ca, a )

a = [ data[3*(i-1)+1:3*i] for i in 1:4]
ca = CellVectorFromDataAndStride(data,3)
test_index_cell_array( ca, a )

p = Point(2.0,1.0)
data = [2,3,1,3,4,4,3,2,5,4,3,4]
ptrs = [1,4,4,7,13]
gid_to_val = [p, 2*p, 3*p, -p, p]
lid_to_gid = CellVectorFromDataAndPtrs(data,ptrs)
ca = CellVectorFromLocalToGlobal(lid_to_gid,gid_to_val)
a = [ gid_to_val[gids] for gids in lid_to_gid ]
test_index_cell_array( ca, a )

b = fill([2,4],(3,4))
c = CartesianArray(b)
lid_to_gid = CellValueFromArray(c)
test_index_cell_array( lid_to_gid, c )
gid_to_val = [1, 2, 3, -1, 2]
ca = CellVectorFromLocalToGlobal(lid_to_gid,gid_to_val)
aa = fill([2,-1],(3,4))
a = CartesianArray(aa)
test_index_cell_array( ca, a )

data = [-2,3,-1,3,4,-4,-3,2,5,4,3,-4]
ptrs = [1,4,4,7,13]
p = 1
gid_to_val_pos = [p, 2*p, 3*p, -p, p]
gid_to_val_neg = [2*p, p, 5*p, 3*p, p]
lid_to_gid = CellVectorFromDataAndPtrs(data,ptrs)
ca = CellVectorFromLocalToGlobalPosAndNeg(
  lid_to_gid,gid_to_val_pos,gid_to_val_neg)
a = [ [1, 3, 2], Int64[], [3, -1, 3], [5, 2, 1, -1, 3, 3] ]
test_index_cell_array( ca, a )

b = fill([2,-4],(3,4))
c = CartesianArray(b)
lid_to_gid = CellValueFromArray(c)
test_index_cell_array( lid_to_gid, c )
gid_to_val_pos = [1, 2, 3, -1, 2]
gid_to_val_neg = [-1, -3, -3, 1, -4]
ca = CellVectorFromLocalToGlobalPosAndNeg(
  lid_to_gid,gid_to_val_pos,gid_to_val_neg)
aa = fill([2,1],(3,4))
a = CartesianArray(aa)
test_index_cell_array( ca, a )

cell_to_x_l = [2,3,1,3,4,4,3,2,5,4,3,4]
cell_to_x_p = [1,4,4,7,13]
cell_to_x = CellVectorFromDataAndPtrs(cell_to_x_l,cell_to_x_p)
x_to_vals_l = [5,4,1,2,3,6,7,8,9,10]
x_to_vals_p = [1,3,5,6,10,11]
x_to_vals = CellVectorFromDataAndPtrs(x_to_vals_l,x_to_vals_p)
cell_to_vals = CellVectorByComposition(cell_to_x, x_to_vals)
a = [
  [1, 2, 3, 5, 4],
  Int[],
  [3, 6, 7, 8, 9, 6, 7, 8, 9],
  [3, 1, 2, 10, 6, 7, 8, 9, 3, 6, 7, 8, 9]]
test_index_cell_array( cell_to_vals, a )

end # module WrappersTests

