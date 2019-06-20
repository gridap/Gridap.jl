module CellValuesAppendTests

using Test
using Gridap
using Gridap.CellValuesGallery

# append

vc1_l = [1, 1, 2, 2, 1, 1, 2, 2]
vc1_p = [1, 2, 4, 4, 5, 6, 8, 9]
vc1 = CellVectorFromDataAndPtrs(vc1_l, vc1_p)

vc2_l = [1, 1, 1, 1, 2, 2, 2, 2]
vc2_p = [1, 2, 2, 4, 7, 8, 9]
vc2 = CellVectorFromDataAndPtrs(vc2_l, vc2_p)

vc3_l = [3, 3, 4, 5, 6]
vc3_p = [1, 2, 6]
vc3 = CellVectorFromDataAndPtrs(vc3_l, vc3_p)

vc = append(vc1,vc2,vc3)

a1 = collect(vc1)
a2 = collect(vc2)
a3 = collect(vc3)
a = vcat(a1,a2,a3)

test_index_cell_array(vc,a)

# local_append

vc1_l = [1, 1, 2, 2, 1, 1, 2, 2]
vc1_p = [1, 2, 4, 4, 5, 6, 8, 9]
vc1 = CellVectorFromDataAndPtrs(vc1_l, vc1_p)

vc2_l = [1, 1, 1, 1, 2, 2, 2, 2]
vc2_p = [1, 2, 2, 4, 6, 7, 8, 9]
vc2 = CellVectorFromDataAndPtrs(vc2_l, vc2_p)

vc3_l = [5, 9, 1, 3, 2, 2, 3, 2]
vc3_p = [1, 1, 2, 4, 5, 7, 8, 9]
vc3 = CellVectorFromDataAndPtrs(vc3_l, vc3_p)

o = 30
vc = local_append(o,vc1,vc2)
a1 = collect(vc1)
a2 = collect(vc2)
a3 = collect(vc3)
a = [ vcat(a1i,a2i.+o) for (a1i,a2i) in zip(a1,a2) ]
test_index_cell_array(vc,a)

o2 = 24 
o3 = 14
vc = local_append((o2,o3),vc1,vc2,vc3)
a = [ vcat(a1i,a2i.+o2,a3i.+o3) for (a1i,a2i,a3i) in zip(a1,a2,a3) ]
test_index_cell_array(vc,a)

end # module
