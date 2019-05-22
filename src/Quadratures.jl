module Quadratures

# Dependencies of this module

using QuadGK: gauss
using Gridap: flatten
using Gridap.Helpers
using Gridap.FieldValues
using Gridap.Polytopes

# Functionality provided by this module

export Quadrature
export TensorProductQuadrature
import Gridap: coordinates, weights

# Abstract types and interfaces

"""
Abstract type representing a quadrature rule on a Polytope in a space of D dimensions
"""
abstract type Quadrature{D} end

coordinates(::Quadrature{D} where D)::Array{Point{D},1} = @abstractmethod

weights(::Quadrature)::Array{Float64,1} = @abstractmethod

# Factories

"""
Factory function to create Quadrature objects in a convenient way
"""
function Quadrature(extrusion::NTuple{D,Int};order::Int) where D
  _quadrature(extrusion,order)
end

# Concrete structs

"""
Tensor product quadrature rule (nodes and weights) on a hyper cube [-1,1]^D
"""
struct TensorProductQuadrature{D} <: Quadrature{D}
  coords::Array{Point{D},1}
  weights::Array{Float64,1}
end

function TensorProductQuadrature(;orders::NTuple{D,Int}) where D
    @assert D == length(orders)
    npoints = [ ceil(Int,(orders[i]+1.0)/2.0) for i in 1:D ]
    quads = [ gauss( eltype(Point{D}), npoints[i] ) for i in 1:D ]
    (coords, weights) = _tensor_product(quads,npoints,Val(D))
    TensorProductQuadrature{D}(coords,weights)
end

coordinates(self::TensorProductQuadrature) = self.coords

weights(self::TensorProductQuadrature) = self.weights

# Helpers

function _quadrature(extrusion::NTuple{D,Int},order) where D
  @notimplementedif any(extrusion == TET_AXIS)
  TensorProductQuadrature(orders=tuple(fill(order,D)...))
end

function _tensor_product(quads,npoints,::Val{D}) where D
  @assert length(quads) == D
  @assert length(npoints) == D
  coords = Array{Point{D},D}(undef,tuple(npoints...))
  weights = Array{Float64,D}(undef,tuple(npoints...))
  _tensor_product!(quads,coords,weights)
  (flatten(coords), flatten(weights))
end

function _tensor_product!(quads,coords::Array{Point{D},D},weights) where D
  p = MPoint{D}(zeros(D))
  for ci in CartesianIndices(coords)
    p[:] = 0.0
    w = 1.0
    for d in 1:D
      xi = quads[d][1][ci[d]]
      wi = quads[d][2][ci[d]]
      p[d] = xi
      w *= wi
    end
    coords[ci] = p
    weights[ci] = w
  end
end

end # module Quadratures
