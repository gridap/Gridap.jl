module QuadratureFactoriesTests

using Test
using Gridap

order = 5

p = Polytope(HEX_AXIS,HEX_AXIS)
quad = Quadrature(p,order)
@test isa(quad,TensorProductQuadrature)


p = Polytope(TET_AXIS,TET_AXIS,TET_AXIS)
quad = Quadrature(p,order)
@test isa(quad,DuffyQuadrature)

end # module
