module MultiCellArraysTests

using Test
using Gridap
using Gridap.CachedArrays

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

mca2 = MultiCellVector([sca,scb],[(1,),(3,)])

sa = [sv sv; sv sv]
sca = ConstantCellArray(sa,l)

sb = [2*sv sv; 3*sv sv]
scb = ConstantCellArray(sb,l)

mca2 = MultiCellMatrix([sca,scb],[(1,5),(3,2)])

end # module MultiCellArraysTests
