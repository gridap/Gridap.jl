module Quadratures

# Dependencies of this module

using QuadGK: gauss
using Gridap
using Gridap.Helpers
using StaticArrays

# Functionality provided by this module

export Quadrature
export TensorProductQuadrature
export coordinates
export weights

# Abstract types and interfaces

"""
Abstract type representing a quadrature rule on a Polytope in a space of D dimensions
"""
abstract type Quadrature{D} end

coordinates(::Quadrature{D} where D)::Array{Point{D,Float64},1} = @abstractmethod

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
  coords::Array{Point{D,Float64},1}
  weights::Array{Float64,1}
end

function TensorProductQuadrature(;orders::NTuple{D,Int}) where D
    @assert D == length(orders)
    npoints = [ ceil(Int,(orders[i]+1.0)/2.0) for i in 1:D ]
    quads = [ gauss( eltype(Point{D,Float64}), npoints[i] ) for i in 1:D ]
    for i in 1:D
      quads[i][1] .+= 1; quads[i][1] .*= 1.0/2.0
      quads[i][2] .*= 1.0/2.0
    end
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
  cis = CartesianIndices(tuple(npoints...))
  n = prod(npoints)
  coords = Vector{Point{D,Float64}}(undef,n)
  weights = Vector{Float64}(undef,n)
  _tensor_product!(quads,coords,weights,cis,D)
  (coords, weights)
end

function _tensor_product!(quads,coords,weights,cis,D)
  p = zero(MVector{D,Float64})
  lis = LinearIndices(cis)
  for ci in cis
    p[:] = 0.0
    w = 1.0
    for d in 1:D
      xi = quads[d][1][ci[d]]
      wi = quads[d][2][ci[d]]
      p[d] = xi
      w *= wi
    end
    li = lis[ci]
    coords[li] = p
    weights[li] = w
  end
end

end # module Quadratures
