module ConstantCellValuesTests

using Test
using LinearAlgebra: inv, det
using Gridap
using Gridap.CellValues
using Gridap.CellValues.ConstantCellValues
using Gridap.CellValues.Testers

l = 10
sv = 1.0
scv = ConstantCellValue(sv,l)

sv2 = 1.1
scv2 = ConstantCellValue(sv2,l)

test_index_cell_value( scv, fill(sv,l) )

sa = [sv, sv, sv]
sca = ConstantCellArray(sa,l)

test_index_cell_array( sca, fill(sa,l) )

sca = ConstantCellVector(sa,l)

test_index_cell_array( sca, fill(sa,l) )

sa2 = [sv sv; sv sv; sv sv]
sca2 = ConstantCellArray(sa2,l)

test_index_cell_array( sca2, fill(sa2,l) )

sca2 = ConstantCellMatrix(sa2,l)

test_index_cell_array( sca2, fill(sa2,l) )

for op in (:+,:-,:*,:/,:(outer),:(inner))
  @eval begin
    scv3 = $op(scv,scv2)
    @test isa(scv3,ConstantCellValue{Float64})
    test_index_cell_value( scv3, fill($op(sv,sv2),l) )
  end
end

scv3 = apply(-,scv,scv2)
@test isa(scv3,ConstantCellValue{Float64})
test_index_cell_value( scv3, fill(-(sv,sv2),l) )

for op in (:+,:-,:*,:/,:(outer),:(inner))
  @eval begin
    sca3 = $op(sca,sca2)
    @test isa(sca3,ConstantCellArray{Float64,2})
    test_index_cell_array( sca3, fill(broadcast($op,sa,sa2),l) )
  end
end

for op in (:+,:-,:*,:/,:(outer),:(inner))
  @eval begin
    sca3 = $op(scv,sca2)
    @test isa(sca3,ConstantCellArray{Float64,2})
    test_index_cell_array( sca3, fill(broadcast($op,sv,sa2),l) )
  end
end

for op in (:+,:-,:(det),:(inv))
  @eval begin
    sca3 = $op(sca2)
    @test isa(sca3,ConstantCellArray{Float64,2})
    test_index_cell_array( sca3, fill(broadcast($op,sa2),l) )
  end
end

sca3 = apply(-,sca2)
@test isa(sca3,ConstantCellArray{Float64,2})
test_index_cell_array( sca3, fill(broadcast(-,sa2),l) )

sca3 = cellsum(sca2,dim=2)
@test isa(sca3,ConstantCellArray{Float64,1})
test_index_cell_array( sca3, fill(reshape(sum(sa2,dims=2),(3,)),l) )

sca3 = cellsum(sca,dim=1)
@test isa(sca3,ConstantCellValue{Float64})
test_index_cell_value( sca3, fill(sum(sa),l) )

sca3 = cellnewaxis(sca2,dim=2)
@test isa(sca3,ConstantCellArray{Float64,3})
test_index_cell_array( sca3, fill(reshape(sa2,(3,1,2)),l) )

end # module ConstantCellValuesTests
