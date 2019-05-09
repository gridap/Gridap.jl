module ConstantCellMapsTests

using Test
using Numa
using Numa.CellValues
using Numa.CellValues.ConstantCellValues
using Numa.Maps
using Numa.CellMaps
using Numa.CellMaps.ConstantCellMaps
using Numa.CellMaps.Testers
using Numa.CellMaps.CellMapValues
using Numa.FieldValues

include("../CellValuesTests/Mocks.jl")
include("../MapsTests/MockMap.jl")
include("../MapsTests/MockBasis.jl")

l = 10
p0 = Point(1.4,2.0)
m = MockMap(p0)
cm = ConstantCellMap(m,l)

p1 = Point(1.0,2.0)
p2 = Point(1.1,2.0)
p3 = Point(1.1,2.3)
ps = [p1,p2,p3]
rs = evaluate(m,ps)
crs = [ rs for i in 1:l]

cv = ConstantCellVector(ps,l)
test_cell_map_without_gradient(cm,cv,crs)
ca = evaluate(cm,cv)
@test isa(ca,ConstantCellArray)

cv = TestCellArray(ps,l)
test_cell_map_without_gradient(cm,cv,crs)
ca = evaluate(cm,cv)
@test isa(ca,CellMapValue)

for op in (:+, :-)
  @eval begin
    ucm = $op(cm)
    test_cell_map_without_gradient(ucm,cv,$op(crs))
    @test isa(ucm,ConstantCellMap)
  end
end

m1 = MockMap(p1)
m2 = MockMap(p2)
cm1 = ConstantCellMap(m1,l)
cm2 = ConstantCellMap(m2,l)
rs1 = evaluate(m1,ps)
rs2 = evaluate(m2,ps)
crs1 = [ rs1 for i in 1:l]
crs2 = [ rs2 for i in 1:l]

for op in (:+, :-, :(inner), :(outer))
  @eval begin
    ucm = $op(cm1,cm2)
    ucrs = [ $op.(rs1,rs2)  for i in 1:l ]
    test_cell_map_without_gradient(ucm,cv,ucrs)
    @test isa(ucm,ConstantCellMap)
  end
end

end # module ConstantCellMapsTests
