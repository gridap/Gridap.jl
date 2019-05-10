module OperationsTests

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

import Numa: gradient

include("../MapsTests/MockMap.jl")
include("../MapsTests/MockBasis.jl")
include("Mocks.jl")

l = 10

p1 = Point(1.0,2.0)
p2 = Point(1.1,2.0)
p3 = Point(1.1,2.3)
ps = [p1,p2,p3]
cv = TestCellArray(ps,l)

p0 = Point(1.4,2.0)
m = MockMap(p0)
@test isa(m,Field)
cm = TestCellMap(m,l)
@test isa(cm,CellField)

rs = evaluate(m,ps)
grs = evaluate(gradient(m),ps)
crs = [ rs for i in 1:l]
gcrs = [ grs for i in 1:l]
test_cell_map_with_gradient(cm,cv,crs,gcrs)

for op in (:+, :-)
  @eval begin
    ucm = $op(cm)
    test_cell_map_without_gradient(ucm,cv,$op(crs))
  end
end

end # module OperationsTests
