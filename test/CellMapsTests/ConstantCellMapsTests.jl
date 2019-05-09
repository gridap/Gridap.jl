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

cv = ConstantCellVector(ps,l)
test_cell_map_without_gradient(cm,cv,[ rs for i in 1:l])
ca = evaluate(cm,cv)
@test isa(ca,ConstantCellArray)

cv = TestCellArray(ps,l)
test_cell_map_without_gradient(cm,cv,[ rs for i in 1:l])
ca = evaluate(cm,cv)
@test isa(ca,CellMapValue)

end # module ConstantCellMapsTests
