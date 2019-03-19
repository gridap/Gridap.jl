using QuadGK

export Quadrature, coordinates, weights
export TensorProductQuadrature, TensorProductQuadratureOld

"""
Abstract type representing a quadrature rule on a Polytope in a space of D dimensions
"""
abstract type Quadrature{D} end

coordinates(::Quadrature{D} where D)::Array{Point{D},1} = @abstractmethod

weights(::Quadrature)::Array{Float64,1} = @abstractmethod

#Concrete implementations

"""
Tensor product quadrature rule (nodes and weights) on a hyper cube [-1,1]^D
"""
struct TensorProductQuadrature{D} <: Quadrature{D}

  coords::Array{Point{D},1}
  weights::Array{Float64,1}

  function TensorProductQuadrature{D}(;orders::Array{Int,1}) where D
    @assert D == length(orders)
    npoints = [ ceil(Int,(orders[i]+1.0)/2.0) for i in 1:D ]
    quads = [ gauss( eltype(Point{D}), npoints[i] ) for i in 1:D ]
    (coords, weights) = tensor_product(quads,npoints,Val(D))
    new{D}(coords,weights)
  end

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


"""
Tensor product quadrature rule (nodes and weights) integrating exactly 2´order´-1 polynomials
"""
struct TensorProductQuadratureOld
  points::Array{Array{Float64,1},1}
  weights::Array{Array{Float64,1},1}
  tppoints::Array{Float64,2}
  tpweights::Array{Float64,1}
  function TensorProductQuadratureOld(gps::Array{Int64,1})
    spdims = length(gps)
    points=Array{Array{Float64,1}}(undef,length(gps))
    weights=Array{Array{Float64,1}}(undef,length(gps))
    points = [gauss(gps[i]) for i=1:length(gps)]
    a = [points[i][1] for i=1:length(gps)]
    b = [points[i][2] for i=1:length(gps)]
    tpa = Array{Float64,spdims+1}(undef,tuple(gps...,spdims))
    tensorfill!(tpa,a)
    tpa = reshape(tpa,:,spdims)
    tpb = Array{Float64,spdims+1}(undef,tuple(gps...,spdims))
    tensorfill!(tpb,b)
    tpb = reshape(tpb,:,spdims)
    tpb = [prod(tpb[i,:]) for i=1:size(tpb,1)]
    return new(a,b,tpa,tpb)
  end
end

@generated function tensorfill!(A::Array{T,N},c) where {T,N}
  quote
    @nloops $N i A begin
      t = @ntuple $N i; d = length(t)
      dim = t[d]
      nodedim = t[dim]
      (@nref $N A i) = c[dim][nodedim]
    end
  end
end

