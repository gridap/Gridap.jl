module MultiFieldCellArraysTests

using Test
using Gridap.MultiField
using Gridap.MultiField: MultiFieldCellArray

l = 10
a = [rand(2,3) for i in 1:l]
b = [rand(2,2) for i in 1:l]
c = [rand(3,3) for i in 1:l]

block_ids = [(1,1),(1,2),(2,1)]

mca = MultiFieldCellArray((a,b,c),block_ids)
mcb = collect(mca)

end # module
