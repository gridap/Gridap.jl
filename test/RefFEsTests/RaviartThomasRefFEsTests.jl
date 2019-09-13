module RaviartThomasRefFEsTests
##

using Gridap
using Gridap.Helpers

using Gridap.RefFEs.RaviartThomasRefFEs

p = Polytope(1,1)

order = 1

reffe = RaviartThomasRefFE(p,order)

dofs = dofbasis(reffe)

p = polytope(reffe)

shb = shfbasis(reffe)

nfdofs = nfacedofs(reffe)

order = 2

reffe = RaviartThomasRefFE(p,order)

dofs = dofbasis(reffe)

p = polytope(reffe)

shb = shfbasis(reffe)

nfdofs = nfacedofs(reffe)

##
end # module
