function TensorProductQuadrature{D}(;orders::Array{Int64,1}) where D
    @assert D == length(orders)
    npoints = [ ceil(Int,(2*orders[i]+1.0)/2.0) for i in 1:D ]
    quads = [ gauss( eltype(Point{D}), npoints[i] ) for i in 1:D ]
    # quads = [ gauss( eltype(Point{D}), orders[i] ) for i in 1:D ]
    # (coords, weights) = tensor_product(quads,npoints,Val(D))
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
