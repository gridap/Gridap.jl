include("MultiCellArrays.jl")

module MultiCellArraysTests

using Test

using ..MultiCellArrays
using Gridap
using Gridap.CellValues
using Gridap.CachedArrays
using Gridap.CellValues.ConstantCellValues

l = 10
sv = 1.0
sa = [sv, sv, sv]
sca = ConstantCellArray(sa,l)

sb = [2*sv, sv, 3*sv]
scb = ConstantCellArray(sb,l)

mca = MultiCellArray([sca,scb],[(1,),(3,)])

@test length(mca) == l

i_to_field = mca.fieldids

for ma in mca
  for (i,a) in enumerate(ma)
    ifield, = i_to_field[i]
    @assert isa(a,CachedArray)
  end
  @assert length(ma) == 2
  @assert isa(ma[1],CachedArray)
  @assert isa(ma[2],CachedArray)
end

end # module MultiCellArraysTests
