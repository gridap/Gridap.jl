module RaviartThomasRefFEsTests
##

using Gridap
using Gridap.Helpers

# using Gridap.RefFEs.RaviartThomasRefFEs

using LinearAlgebra

using Test

using Gridap.RefFEs.RTRefFEs
#####
# Create a general Ref FE constructor, which generalizes the current implementations
# for nodal (Lagrangian) and non-nodal (RaviartThomas) RefFEs.
# * A method that generates the geomap from ref n-face to n-face in polytope
# * A functor to be provided by the user (optionally) for every n-face dimension
#####
##
p = Polytope(1,1)
order = 3
ref_fe = RTRefFE(p,3)
# ref_fe2 = RTRefFE(p,3)
# ref_fe.shfbasis.changeofbasis == ref_fe2.shfbasis.changeofbasis

# Cell moments


for order in 1:3

  p = Polytope(1,1)

  reffe = RTRefFE(p,order)

  dofsb = dofbasis(reffe)

  p = polytope(reffe)

  shb = shfbasis(reffe)

  nfdofs = nfacedofs(reffe)


  b = dofsb
  shfs = reffe.shfbasis
  kk = evaluate(b,shfs)
  @test kk ≈ Matrix(1.0I, size(kk))

  fun(x) = VectorValue(x[1],2*x[1])
  D = 2
  T = VectorValue{2,Float64}
  field = AnalyticalField(fun,D)

  isa(dofsb,DOFBasis{2,T})
  isa(field,Field{2,T})

  evaluate(b,field)
end
##
end # module
