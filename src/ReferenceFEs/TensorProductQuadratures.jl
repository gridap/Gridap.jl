
"""
    struct TensorProduct <: QuadratureName

Tensor product quadrature rule for n-cubes, obtained as the
tensor product of 1d Gauss-Legendre quadratures.

# Constructor:

    Quadrature(p::Polytope{D}, tensor_product, degrees::Integer; T::Type=Float64)
    Quadrature(p::Polytope{D}, tensor_product, degrees::NTuple{D,Integer}; T::Type=Float64)
"""
struct TensorProduct <: QuadratureName end

"""
    const tensor_product = TensorProduct()
"""
const tensor_product = TensorProduct()

function Quadrature(p::Polytope,::TensorProduct,degrees; T::Type=Float64)
  @assert is_n_cube(p) "Tensor product quadrature rule only for n-cubes."
  D = num_dims(p)
  if iszero(D)
    coords, weights = [zero(Point{0,T})], [one(T)]
  else
    @assert length(degrees) == D
    quad_1d = [ gauss_legendre_quadrature(degree;T) for degree in degrees ]
    coords_1d, weights_1d = map(first, quad_1d), map(last, quad_1d)
    coords, weights = _tensor_product(coords_1d, weights_1d)
  end
  GenericQuadrature(coords,weights,"Tensor product of 1d Gauss-Legendre quadratures of degrees $degrees")
end

function Quadrature(p::Polytope,name::TensorProduct,degree::Integer;kwargs...)
  degrees = ntuple(i->degree,Val(num_dims(p)))
  Quadrature(p,name,degrees;kwargs...)
end

function Quadrature(p::Polytope,::TensorProduct,quadratures::Vector{<:Quadrature{1}})
  @assert is_n_cube(p) "Tensor product quadrature rule only for n-cubes."
  D = num_dims(p)
  @assert length(quadratures) == D

  coords_1d = map(q -> map(xi -> xi[1], get_coordinates(q)), quadratures)
  weights_1d = map(get_weights, quadratures)
  coords, weights = _tensor_product(coords_1d, weights_1d)

  names_1d = join(map(get_name, quadratures)," \n - ")
  GenericQuadrature(coords,weights,"Tensor product of 1d quadratures given by: \n "*names_1d)
end

# Helpers

function _tensor_product(
  d_to_xs::Vector{Vector{T}},
  d_to_ws::Vector{Vector{W}}
) where {T,W}
  D = length(d_to_xs)
  d_to_n = map(length, d_to_xs)
  cis = CartesianIndices(tuple(d_to_n...))

  n = prod(d_to_n)
  xs = Vector{Point{D,T}}(undef,n)
  ws = Vector{W}(undef,n)
  _tensor_product!(xs,ws,d_to_xs,d_to_ws,cis)
  return xs, ws
end

function _tensor_product!(xs,ws,d_to_xs,d_to_ws,cis)
  P = eltype(xs)
  T = eltype(P)
  W = eltype(ws)
  D = length(P)

  k = 1
  z = zero(T)
  p = zeros(T,D)
  for ci in cis
    fill!(p,z)
    w = one(W)
    for d in 1:D
      p[d] = d_to_xs[d][ci[d]]
      w *= d_to_ws[d][ci[d]]
    end
    xs[k] = p
    ws[k] = w
    k += 1
  end
  return xs, ws
end
