
"""
    struct TensorProductQuadrature{D,T} <: Quadrature{D,T}
      coordinates::Vector{Point{D,T}}
      weights::Vector{T}
    end

Tensor product quadrature rule (nodes and weights) on a hyper cube [0,1]^D
"""
struct TensorProductQuadrature{D,T} <: Quadrature{D,T}
  coordinates::Vector{Point{D,T}}
  weights::Vector{T}
end

get_coordinates(q::TensorProductQuadrature) = q.coordinates

get_weights(q::TensorProductQuadrature) = q.weights

"""
    TensorProductQuadrature(degrees::NTuple{D}) where D
    TensorProductQuadrature(degrees::Point{D}) where D
"""
function TensorProductQuadrature(degrees::NTuple{D,Int}) where D
  TensorProductQuadrature{D}(degrees)
end

function TensorProductQuadrature(degrees::Point{D}) where D
  TensorProductQuadrature{D}(degrees)
end

"""
    TensorProductQuadrature{D}(degree::Integer) where D
    TensorProductQuadrature{D}(degrees) where D
"""
function TensorProductQuadrature{D}(degrees) where D
    @assert D == length(degrees)
    T = Float64
    npoints = [ ceil(Int,(degrees[i]+1.0)/2.0) for i in 1:D ]
    quads = [ gauss( eltype(Point{D,T}), npoints[i] ) for i in 1:D ]
    for i in 1:D
      quads[i][1] .+= 1;
      quads[i][1] .*= 1.0/2.0
      quads[i][2] .*= 1.0/2.0
    end
    (coords, weights) = _tensor_product(Point{D,T},quads,npoints)
    TensorProductQuadrature(coords,weights)
end

function TensorProductQuadrature{D}(degree::Integer) where D
  degrees = tfill(degree,Val{D}())
  TensorProductQuadrature(degrees)
end

# Helpers

function _tensor_product(::Type{Point{D,T}},quads,npoints) where {D,T}
  @assert length(quads) == D
  @assert length(npoints) == D
  cis = CartesianIndices(tuple(npoints...))
  n = prod(npoints)
  coords = Vector{Point{D,T}}(undef,n)
  weights = Vector{T}(undef,n)
  _tensor_product!(quads,coords,weights,cis)
  (coords, weights)
end

function _tensor_product(::Type{Point{0,T}},quads,npoints) where T
  coords = [zero(Point{0,T})]
  weights = [one(T)]
  coords, weights
end

function _tensor_product!(quads,coords,weights,cis)
  p = zero(Mutable(eltype(coords)))
  D = length(p)
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
