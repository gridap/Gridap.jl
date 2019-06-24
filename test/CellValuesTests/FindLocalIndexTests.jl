module FindLocalIndexTests

using Gridap
using Gridap.CellValuesGallery

c = [[2,3,4,8],[1,8,4,4,30],[3,5]]
cv = CellValueFromArray(c)

g = [4,30,3]
cn = CellValueFromArray(g)

r = [3,5,1]

li = find_local_index(cn,cv)
test_index_cell_value(li,r)

end # module
