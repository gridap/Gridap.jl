module QuadratureFactories

using Gridap
using Gridap.Helpers

import Gridap: Quadrature

"""
Factory function to create Quadrature objects in a convenient way
"""
function Quadrature(p::Polytope,order::Int)
  Quadrature(p.extrusion.array.data,order=order)
end

function Quadrature(extrusion::NTuple{D,Int};order::Int) where D

  if all(extrusion .== HEX_AXIS)
    return TensorProductQuadrature(orders=tuple(fill(order,D)...))

  elseif all(extrusion .== TET_AXIS)
    return DuffyQuadrature{D}(order)

  else
    @notimplemented
  end

end

end # module
