module CellNumbersTests

using Test
using Gridap

using ..CellValuesMocks

l = 10
a = 3.5

cn = TestIterCellValue(a,l)
@test isa(cn,CellNumber)
@test isa(cn,IterCellNumber)
r = [a for i in 1:l]
test_iter_cell_number(cn,r)
@test !isa(cn,CellArray)
@test !isa(cn,IterCellArray)
@test !isa(cn,CellMap)
@test !isa(cn,IterCellMap)

cn = TestIndexCellValue(a,l)
@test isa(cn,CellNumber)
@test isa(cn,IndexCellNumber)
test_index_cell_number(cn,r)
@test !isa(cn,CellArray)
@test !isa(cn,IndexCellArray)
@test !isa(cn,CellMap)
@test !isa(cn,IndexCellMap)

end # module
