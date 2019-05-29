include("MultiCellArrays.jl")

module MultiCellArraysTests

using Test

using ..MultiCellArrays
using Gridap
using Gridap.CellValues
using Gridap.CachedArrays
using Gridap.CellValues.ConstantCellValues
using Gridap.CellValues.Operations: CellArrayFromBroadcastUnaryOp

l = 10
sv = 1.0
sa = [sv, sv, sv]
sca = ConstantCellArray(sa,l)
sca = CellArrayFromBroadcastUnaryOp(+,sca)

sb = [2*sv, sv, 3*sv]
scb = ConstantCellArray(sb,l)
scb = CellArrayFromBroadcastUnaryOp(+,scb)

T = Float64
N = 1
@test isa(sca,CellValue{CachedArray{T,N,Array{T,N}}})
@test isa(scb,CellValue{CachedArray{T,N,Array{T,N}}})

mca = MultiCellArray([sca,scb],[(1,),(3,)])

@test length(mca) == l

for ma in mca
  @assert length(ma) == 2
  for (a, (ifield,)) in eachblock(ma)
    @show a
    @show ifield
  end
  @show ma[1]
  @show ma[2]
end

end # module MultiCellArraysTests
