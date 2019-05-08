module OperationsTests

using Test
using Numa
using Numa.FieldValues
using Numa.CellValues
using Numa.FieldValues
using LinearAlgebra: det, inv
using StaticArrays
using Numa.CellValues: _custom_broadcast!

include("Helpers.jl")
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

using Numa.CellValues: CellValueFromUnaryOp
using Numa.CellValues: CellValueFromBinaryOp
using Numa.CellValues: CellArrayFromBroadcastUnaryOp
using Numa.CellValues: CellArrayFromCellSum
using Numa.CellValues: CellArrayFromCellNewAxis
using Numa.CellValues: CellValueFromCellArrayReduce
using Numa.CellValues: CellArrayFromBroadcastBinaryOp

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

function _custom_broadcast(op,a,b)
  broadcast(op,a,b)
end

function _custom_broadcast(
  op, a::SArray{Size,T}, b::AbstractArray{S,M}) where {Size,T,S,M}
  R = Base._return_type(op,Tuple{typeof(a),S})
  v = Array{R,M}(undef,size(b))
  _custom_broadcast!(op,v,a,b)
  v
end

function _custom_broadcast(
  op, a::AbstractArray{S,M}, b::SArray{Size,T}) where {Size,T,S,M}
  R = Base._return_type(op,Tuple{S,typeof(b)})
  v = Array{R,M}(undef,size(a))
  _custom_broadcast!(op,v,a,b)
  v
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

end # module OperationsTests

#
#@testset "Operations" begin
#
#
#
#
#  for op in (:+,:-,:*,:/,:(outer),:(inner))
#    @eval begin
#      sca3 = $op(sca,sca2)
#      @test isa(sca3,CellArrayFromBroadcastBinaryOp{typeof($op),Float64,2})
#      @test length(sca3) == l
#      @test cellsize(sca3) == (3,2)
#      for vi in sca3
#        @assert vi == $op.(sa,sa2)
#      end
#    end
#  end
#
#  for op in (:+,:-,:*,:/,:(outer),:(inner))
#    @eval begin
#      sca3 = $op(scv,sca2)
#      @test isa(sca3,CellArrayFromBroadcastBinaryOp{typeof($op),Float64,2})
#      @test length(sca3) == l
#      @test cellsize(sca3) == (3,2)
#      for vi in sca3
#        @assert vi == $op.(sv,sa2)
#      end
#    end
#  end
#
#  for op in (:+,:-,:*,:/,:(outer),:(inner))
#    @eval begin
#      sca3 = $op(sca2,scv)
#      @test isa(sca3,CellArrayFromBroadcastBinaryOp{typeof($op),Float64,2})
#      @test length(sca3) == l
#      @test cellsize(sca3) == (3,2)
#      for vi in sca3
#        @assert vi == $op.(sa2,sv)
#      end
#    end
#  end
#
#  sca3 = cellsum(sca2,dim=2)
#  @test isa(sca3,CellArrayFromCellSum{2,1,typeof(sca2),Float64})
#  @test length(sca3) == l
#  @test cellsize(sca3) == (3,)
#  for vi in sca3
#    @assert vi == reshape(sum(sa2,dims=2),(3,))
#  end
#
#  sca3 = cellsum(sca,dim=1)
#  @test isa(sca3,CellValueFromCellArrayReduce{Float64,typeof(sum),typeof(sca)})
#  @test length(sca3) == l
#  for vi in sca3
#    @assert vi == sum(sa)
#  end
#
#  sca3 = cellmean(sca)
#  @test isa(sca3,CellValueFromCellArrayReduce{Float64,typeof(Numa.CellValues.mean),typeof(sca)})
#  @test length(sca3) == l
#  for vi in sca3
#    @assert vi == Numa.CellValues.mean(sa)
#  end
#
#  sca3 = cellnewaxis(sca2,dim=2)
#  @test isa(sca3,CellArrayFromCellNewAxis{2,typeof(sca2),Float64,3})
#  @test length(sca3) == l
#  @test cellsize(sca3) == (3,1,2)
#  for vi in sca3
#    @assert vi == reshape(sa2,(3,1,2))
#  end
#
#end
#
#@testset "FlattedCellArray" begin
#
#  scv3 = flatten(sca)
#
#  c = collect(scv3)
#
#  i = 1
#  for a in sca
#    for ai in a
#      @assert c[i] == ai
#      i += 1
#    end
#  end
#
#end

