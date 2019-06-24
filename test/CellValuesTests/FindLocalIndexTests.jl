module FindLocalIndexTests

using Gridap
using Gridap.CellValuesGallery

c = [[2,1,3,4,8],[1,8,4,3,4,30],[3,5,2]]
cv = CellValueFromArray(c)

g = [1,3,2,1]
cn = CellValueFromArray(g)

r = [2,3,4,4]

li = find_local_index(cn,cv)
test_index_cell_value(li,r)

li = get_local_item(cv,2)
r = [1,8,5]
test_index_cell_value(li,r)

end # module
