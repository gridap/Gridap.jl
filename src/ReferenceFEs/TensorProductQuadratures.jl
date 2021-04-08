
struct TensorProduct <: QuadratureName end

const tensor_product = TensorProduct()

function Quadrature(p::Polytope,::TensorProduct,degrees)
  @assert is_n_cube(p) """\n
  Tensor product quadrature rule only for n-cubes.
  """
  @assert length(degrees) == num_dims(p)
  _tensor_product_legendre(degrees)
end

function Quadrature(p::Polytope,name::TensorProduct,degree::Integer)
  degrees = ntuple(i->degree,Val(num_dims(p)))
  Quadrature(p,name,degrees)
end

# Low level constructor

function _tensor_product_legendre(degrees)
    D = length(degrees)
    T = Float64
    npoints = [ ceil(Int,(degrees[i]+1.0)/2.0) for i in 1:D ]
    quads = [ gauss( eltype(Point{D,T}), npoints[i] ) for i in 1:D ]
    for i in 1:D
      quads[i][1] .+= 1;
      quads[i][1] .*= 1.0/2.0
      quads[i][2] .*= 1.0/2.0
    end
    (coords, weights) = _tensor_product(Point{D,T},quads,npoints)
    GenericQuadrature(coords,weights,"Tensor product of 1d Gauss-Legendre quadratures of degrees $degrees")
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
