module IdentityCellNumbersTests

using Test
using Gridap
using Gridap.CellValuesGallery

a = [2,3,4,6,1,8,3,5]
cn = CellValueFromArray(a)

l = length(a)
id = IdentityCellNumber(Int,l)
r = collect(1:l)
test_index_cell_number(id,r)

cn2 = reindex(cn,id)

@test cn2 === cn

end # module
