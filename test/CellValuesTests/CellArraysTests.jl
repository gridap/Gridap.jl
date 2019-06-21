module CellArraysTests

using Test
using Gridap

using ..CellValuesMocks

l = 10
a = rand(3,2,4,2)

cn = TestIterCellValue(a,l)
@test isa(cn,CellArray)
@test isa(cn,IterCellArray)
r = [a for i in 1:l]
test_iter_cell_array(cn,r)
@test !isa(cn,CellNumber)
@test !isa(cn,IterCellNumber)
@test !isa(cn,CellMap)
@test !isa(cn,IterCellMap)

cn = TestIndexCellValue(a,l)
@test isa(cn,CellArray)
@test isa(cn,IndexCellArray)
test_index_cell_array(cn,r)
@test !isa(cn,CellNumber)
@test !isa(cn,IndexCellNumber)
@test !isa(cn,CellMap)
@test !isa(cn,IndexCellMap)

a = rand(3)

cn = TestIterCellValue(a,l)
@test isa(cn,CellVector)
@test isa(cn,IterCellVector)
@test !isa(cn,CellNumber)
@test !isa(cn,IterCellNumber)
@test !isa(cn,CellMap)
@test !isa(cn,IterCellMap)

cn = TestIndexCellValue(a,l)
@test isa(cn,CellVector)
@test isa(cn,IndexCellVector)
@test !isa(cn,CellNumber)
@test !isa(cn,IndexCellNumber)
@test !isa(cn,CellMap)
@test !isa(cn,IndexCellMap)

a = rand(3,5)

cn = TestIterCellValue(a,l)
@test isa(cn,CellMatrix)
@test isa(cn,IterCellMatrix)
@test !isa(cn,CellNumber)
@test !isa(cn,IterCellNumber)
@test !isa(cn,CellMap)
@test !isa(cn,IterCellMap)

cn = TestIndexCellValue(a,l)
@test isa(cn,CellMatrix)
@test isa(cn,IndexCellMatrix)
@test !isa(cn,CellNumber)
@test !isa(cn,IndexCellNumber)
@test !isa(cn,CellMap)
@test !isa(cn,IndexCellMap)

end # module
