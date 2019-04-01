module Quadratures

export Quadrature
export TensorProductQuadrature
export coordinates
export weights

using QuadGK: gauss

using Numa # @santiagobadia : To be eliminated
using Numa.Helpers
using Numa.FieldValues

# Abstract types and interfaces

"""
Abstract type representing a quadrature rule on a Polytope in a space of D dimensions
"""
abstract type Quadrature{D} end

coordinates(::Quadrature{D} where D)::Array{Point{D},1} = @abstractmethod

weights(::Quadrature)::Array{Float64,1} = @abstractmethod

#Concrete structs

"""
Tensor product quadrature rule (nodes and weights) on a hyper cube [-1,1]^D
"""
struct TensorProductQuadrature{D} <: Quadrature{D}
  coords::Array{Point{D},1}
  weights::Array{Float64,1}
end

# Methods

function TensorProductQuadrature(;orders::NTuple{D,Int}) where D
    @assert D == length(orders)
    npoints = [ ceil(Int,(orders[i]+1.0)/2.0) for i in 1:D ]
    quads = [ gauss( eltype(Point{D}), npoints[i] ) for i in 1:D ]
    (coords, weights) = tensor_product(quads,npoints,Val(D))
    TensorProductQuadrature{D}(coords,weights)
  end

function tensor_product(quads,npoints,::Val{D}) where D
  @assert length(quads) == D
  @assert length(npoints) == D
  coords = Array{Point{D},D}(undef,tuple(npoints...))
  weights = Array{Float64,D}(undef,tuple(npoints...))
  tensor_product!(quads,coords,weights)
  (flatten(coords), flatten(weights))
end

function tensor_product!(quads,coords::Array{Point{D},D},weights) where D
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

coordinates(self::TensorProductQuadrature) = self.coords

weights(self::TensorProductQuadrature) = self.weights

end # module Quadratures
