module RaviartThomasRefFEsTests
##

using Gridap
using Gridap.Helpers

using Gridap.RefFEs.RaviartThomasRefFEs

using LinearAlgebra

using Test


# p = Polytope(1,1)
# fns, o = facet_normals(p)
# fns
# o
#
# fns[1]*o[1]

for order in 1:3

  p = Polytope(1,1)

  reffe = RaviartThomasRefFE(p,order)

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
