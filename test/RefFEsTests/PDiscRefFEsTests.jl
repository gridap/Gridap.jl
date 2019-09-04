module PDiscRefFEsTests

using Test
using Gridap

D = 2
order = 1

reffe = PDiscRefFE(Float64,D,order)

@test reffe.nfacedofs[end] == [1,2,3]

dofbasis(reffe)
polytope(reffe)
shfbasis(reffe)
nfacedofs(reffe)

end # module
