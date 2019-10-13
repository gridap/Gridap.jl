module GenericRefFEsTests
##

using Gridap
using Gridap.Helpers

# using Gridap.RefFEs.RaviartThomasRefFEs

using LinearAlgebra

using Test

using Gridap.GenericRefFEs

#####
# Create a general Ref FE constructor, which generalizes the current implementations
# for nodal (Lagrangian) and non-nodal (RaviartThomas) RefFEs.
# * A method that generates the geomap from ref n-face to n-face in polytope
# * A functor to be provided by the user (optionally) for every n-face dimension
#####
F = Gridap.RefFEs.GenericRefFEs

p = Polytope(1,1,1)

order = 4

ref2 = NedelecRefFE(p,3)
ref2 = F.MyNedelecRefFE(p,3)
# ref1.shfbasis.changeofbasis == ref2.shfbasis.changeofbasis



# 1. Prebasis
prebasis = CurlGradMonomialBasis(VectorValue{dim(p),Float64},order)

# Nface nodes, moments, and prebasis evaluated at nodes
nf_nodes, nf_moments, pb_moments = _initialize_arrays(prebasis,p)

ccips, cmoments = G._RT_cell_values(p,order)

pbasis_ccips = [evaluate(prebasis,ps) for ps in ccips]

cms_preb = [bps*ms' for (bps,ms) in zip(pbasis_ccips,cmoments)]
F._nfaces_array_dim!(p,dim(p),nf_moments,cmoments)
F._nfaces_array_dim!(p,dim(p),nf_nodes,ccips)
pb_moments = hcat(pb_moments,cms_preb)
pb_moments
cms_preb

############################

function test_reffe(reffe,order)

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
p = Polytope(1,1)

for order in 1:4
  reffe = NedelecRefFE(p,order)
  test_reffe(reffe,order)
end

for order in 1:4
  reffe = RTRefFE(p,order)
  test_reffe(reffe,order)
end

##
end # module
