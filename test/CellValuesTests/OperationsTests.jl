module OperationsTests

using Test
using Gridap
using Gridap.FieldValues
using Gridap.CellValues
using Gridap.FieldValues
using LinearAlgebra: det, inv
using StaticArrays
using Gridap.CellValues.Operations: _custom_broadcast
using Gridap.CellValues.Operations: mean
using Gridap.CellValues.Testers

include("Mocks.jl")

l = 10
sv = 1.0
scv = TestCellValue(sv,l)
test_index_cell_value(scv,fill(sv,l))

sv_t = TensorValue(1.0,2.0,0.0,1.0)
scv_t = TestCellValue(sv_t,l)
test_index_cell_value(scv_t,fill(sv_t,l))

sv_v = VectorValue(1.1,2.3)
scv_v = TestCellValue(sv_v,l)

sv2 = 1.1
scv2 = TestCellValue(sv2,l)

sv2_t = TensorValue(0.1,2.0,1.0,1.0)
scv2_t = TestCellValue(sv2_t,l)

sv2_v = VectorValue(1.1,2.0)
scv2_v = TestCellValue(sv2_v,l)

sa = [sv, sv, sv]
sca = TestCellArray(sa,l)
test_index_cell_array(sca,fill(sa,l))

sa2 = [sv sv; sv sv; sv sv]
sca2 = TestCellArray(sa2,l)

sa2_v = fill(sv_v,(3,2))
sca2_v = TestCellArray(sa2_v,l)
test_index_cell_array(sca2_v,fill(sa2_v,l))

sa_v = [sv_v, sv_v, sv_v]
sca_v = TestCellArray(sa_v,l)
test_index_cell_array(sca_v,fill(sa_v,l))

sa_t = [sv_t, sv_t, sv_t]
sca_t = TestCellArray(sa_t,l)
test_index_cell_array(sca_t,fill(sa_t,l))

sa2_t = fill(sv_t,(3,2))
sca2_t = TestCellArray(sa2_t,l)
test_index_cell_array(sca2_t,fill(sa2_t,l))

@test eltype(scv) == Float64

s = string(sca)
s0 = """
1 -> [1.0, 1.0, 1.0]
2 -> [1.0, 1.0, 1.0]
3 -> [1.0, 1.0, 1.0]
4 -> [1.0, 1.0, 1.0]
5 -> [1.0, 1.0, 1.0]
6 -> [1.0, 1.0, 1.0]
7 -> [1.0, 1.0, 1.0]
8 -> [1.0, 1.0, 1.0]
9 -> [1.0, 1.0, 1.0]
10 -> [1.0, 1.0, 1.0]
"""
@test s == s0

@test cellsize(sca2) == (3,2)
@test cellsize(sca2,1) == 3
@test cellsize(sca2,2) == 2
@test celllength(sca2) == 6

using Gridap.CellValues.Operations: CellValueFromUnaryOp
using Gridap.CellValues.Operations: CellValueFromBinaryOp
using Gridap.CellValues.Operations: CellArrayFromBroadcastUnaryOp
using Gridap.CellValues.Operations: CellArrayFromCellSum
using Gridap.CellValues.Operations: CellArrayFromCellNewAxis
using Gridap.CellValues.Operations: CellValueFromCellArrayReduce
using Gridap.CellValues.Operations: CellArrayFromBroadcastBinaryOp

for op in (:+,:-,:(inv),:(det))
  @eval begin
    for scvi in [ scv, scv_t ]
      svi = scvi.a
      scv3 = $op(scvi)
      isa(scv3,CellValueFromUnaryOp)
      test_iter_cell_value(scv3, fill($op(svi),l))
    end
  end
end

for op in (:+,:-,:*,:/,:(inner))
  @eval begin
    for (scvi,scv2i) in [ (scv,scv2) , (scv_t, scv2_t) ]
      svi = scvi.a
      svi2 = scv2i.a
      scv3 = $op(scvi,scv2i)
      @test isa(scv3,CellValueFromBinaryOp)
      test_iter_cell_value(scv3, fill($op(svi,svi2),l))
    end
  end
end

scv3 = apply(+,scv,scv2)
test_iter_cell_value(scv3, fill(+(scv.a,scv2.a),l))

for op in (:+,:-,:(inv),:(det))
  @eval begin
    for scai in [ sca, sca_t ]
      sai = scai.a
      sca3 = $op(scai)
      @test isa(sca3,CellArrayFromBroadcastUnaryOp)
      test_iter_cell_array(sca3,fill($op.(sai),l))
    end
  end
end

for op in (:+,:-,:*,:/,:(outer),:(inner))
  @eval begin
    for (i,scai) in enumerate([ sca, scv ])
      sai = scai.a
      for (j,scaj) in enumerate([sca2, scv2])
        saj = scaj.a
        sca3 = $op(scai,scaj)
        if (i,j) != (2,2)
          @test isa(sca3,CellArrayFromBroadcastBinaryOp)
          test_iter_cell_array(sca3,fill($op.(sai,saj),l))
        else
          @test isa(sca3,CellValueFromBinaryOp)
          test_iter_cell_value(sca3,fill($op.(sai,saj),l))
        end
      end
    end
  end
end

for op in (:+,:-,:*,:/,:(inner))
  @eval begin
    for (i,scai) in enumerate([ sca_t, scv_t])
      sai = scai.a
      for (j,scaj) in enumerate([sca2_t, scv2_t])
        saj = scaj.a
        sca3 = $op(scai,scaj)
        if (i,j) != (2,2)
          @test isa(sca3,CellArrayFromBroadcastBinaryOp)
          test_iter_cell_array(sca3,fill(_custom_broadcast($op,sai,saj),l))
        else
          @test isa(sca3,CellValueFromBinaryOp)
          test_iter_cell_value(sca3,fill($op(sai,saj),l))
        end
      end
    end
  end
end

for scai in [sca2, sca2_v, sca2_t]
  sai = scai.a
  sca3 = cellsum(scai,dim=2)
  @test isa(sca3,CellArrayFromCellSum{2})
  test_iter_cell_array(sca3,fill(reshape(sum(sai,dims=2),(3,)),l))
end

for scai in [sca, sca_v, sca_t]
  sai = scai.a
  sca3 = cellsum(scai,dim=1)
  @test isa(sca3,CellValueFromCellArrayReduce)
  test_iter_cell_value(sca3,fill(sum(sai),l))
end

for scai in [sca, sca_v, sca_t]
  sai = scai.a
  sca3 = cellmean(scai)
  @test isa(sca3,CellValueFromCellArrayReduce)
  test_iter_cell_value(sca3,fill(mean(sai),l))
end

for scai in [sca2, sca2_v, sca2_t]
  sai = scai.a
  sca3 = cellnewaxis(scai,dim=2)
  @test isa(sca3,CellArrayFromCellNewAxis{2})
  test_iter_cell_array(sca3,fill(reshape(sai,(3,1,2)),l))
end

end # module OperationsTests

