module CompressedCellValuesTests

using Test
using Gridap
using Gridap.CellValuesGallery
using ..MapsMocks

values = [10,20,31]
ptrs = [1,2,3,3,2,2]
r = values[ptrs]
cv = IterCompressedCellValue(values,ptrs)
test_iter_cell_value(cv,r)

cv = IndexCompressedCellValue(values,ptrs)
test_index_cell_value(cv,r)

l = 10
v = 3
ccv = ConstantCellValue(v,l)
cv = IndexCompressedCellValue(ccv)
test_index_cell_value(cv,fill(v,l))

@test cv == cv
@test cv ≈ cv

values = [10,20,31]
ptrs = [1,2,3,3,2,2]
cv = CompressedCellValue(values,ptrs)

cv2 = apply(-,cv)
@test isa(cv2,IndexCompressedCellValue)
r = -values[ptrs]
test_index_cell_value(cv2,r)

values1 = [10,20,31]
values2 = [11,21,51]
ptrs = [1,2,3,3,2,2]
cv1 = CompressedCellValue(values1,ptrs)
cv2 = CompressedCellValue(values2,ptrs)
cv3 = apply(-,cv1,cv2)
@test isa(cv3,IndexCompressedCellValue)
r = values1[ptrs] - values2[ptrs]
test_index_cell_value(cv3,r)

values1 = [10,20,31]
values2 = [11,21,51]
ptrs1 = [1,2,3,3,2,2]
ptrs2 = [1,1,3,2,3,2]
cv1 = CompressedCellValue(values1,ptrs1)
cv2 = CompressedCellValue(values2,ptrs2)
cv3 = apply(-,cv1,cv2)
@test ! isa(cv3,IndexCompressedCellValue)
r = values1[ptrs1] - values2[ptrs2]
test_index_cell_value(cv3,r)

values1 = [10,20,31]
values2 = [11,21,51]
ptrs1 = [1,2,3,3,2,2]
ptrs2 = [1,1,3,2,3,2]
cv1 = IterCompressedCellValue(values1,ptrs1)
cv2 = IterCompressedCellValue(values2,ptrs2)
cv3 = apply(-,cv1,cv2)
@test ! isa(cv3,IndexCompressedCellValue)
r = values1[ptrs1] - values2[ptrs2]
test_iter_cell_value(cv3,r)

values = [[10,30],[20,10,40],[31,]]
ptrs = [1,2,3,3,2,2]
cv = CompressedCellValue(values,ptrs)

cv2 = apply(-,cv,broadcast=true)
@test isa(cv2,IndexCompressedCellValue)
r = -values[ptrs]
test_index_cell_array(cv2,r)

values1 = [[10,33],[20,12,40],[71,]]
values2 = [[10,30],[20,10,40],[31,]]
ptrs = [1,2,3,3,2,2]
cv1 = CompressedCellValue(values1,ptrs)
cv2 = CompressedCellValue(values2,ptrs)
cv3 = apply(-,cv1,cv2,broadcast=true)
@test isa(cv3,IndexCompressedCellValue)
r = values1[ptrs] - values2[ptrs]
test_index_cell_array(cv3,r)

values1 = [[10,33,2],[20,12,40],[71,1,4]]
values2 = [[10,30,3],[20,10,40],[31,3,5]]
ptrs1 = [1,3,2,3,2,2]
ptrs2 = [1,2,3,3,2,2]
cv1 = CompressedCellValue(values1,ptrs1)
cv2 = CompressedCellValue(values2,ptrs2)
cv3 = apply(-,cv1,cv2,broadcast=true)
r = values1[ptrs1] - values2[ptrs2]
@test ! isa(cv3,IndexCompressedCellValue)
test_index_cell_array(cv3,r)


a = VectorValue(10,10)
b = VectorValue(15,20)
p1 = VectorValue(1,1)
p2 = VectorValue(2,2)
p3 = VectorValue(3,3)
p = [p1,p2,p3]
m = MockMap(a)

l = 5
cm = ConstantCellMap(m,l)
vals = [[p1,p3], [p1,p3,p2], ]
ptrs = [1,2,1,1,2]
ca = CompressedCellValue(vals,ptrs)
cv = evaluate(cm,ca)
@test isa(cv,CompressedCellValue)
rvals = [ evaluate(m,pj) for pj in vals  ]
cr = rvals[ptrs]
test_index_cell_map_with_index_arg(cm,ca,cr)

m1 = MockMap(a)
m2 = MockMap(b)
mvals = [m1,m2]
pvals = [[p1,p3], [p1,p3,p2], ]
ptrs = [1,2,1,1,2]
cm = CompressedCellMap(mvals,ptrs)
ca = CompressedCellArray(pvals,ptrs)
cv = evaluate(cm,ca)
@test isa(cv,CompressedCellValue)
rvals = [ evaluate(mi,ai) for (mi,ai) in zip(mvals,pvals) ]
cr = CompressedCellValue(rvals,ptrs)
test_index_cell_map_with_index_arg(cm,ca,cr)

values = [10,20,31]
ptrs = [1,2,3,3,2,2]
inds = [2,1,4]
cv = CompressedCellValue(values,ptrs)
indices = CellValueFromArray(inds)
cv2 = reindex(cv,indices)
@test isa(cv2,CompressedCellValue)
r = values[ptrs[inds]]
test_index_cell_value(cv2,r)

end # module
