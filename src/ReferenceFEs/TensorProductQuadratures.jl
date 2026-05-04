
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

function Quadrature(
  p::Polytope, ::TensorProduct, degrees; T::Type{<:AbstractFloat}=Float64
)
  msg = """\n
  `tensor_product` quadrature rule only for n-cubes.
  """
  @assert is_n_cube(p) msg

  msg = """\n
  Incompatible degrees. Please specify one degree per dimension.
  """
  @assert (length(degrees) == num_dims(p))

  coords, weights = _tensor_product_legendre(degrees; T)
  name = "Tensor product of 1d Gauss-Legendre quadratures of degrees $degrees"
  GenericQuadrature(coords, weights, name)
end

function Quadrature(
  p::Polytope, name::TensorProduct, degree::Integer; T::Type{<:AbstractFloat}=Float64
)
  degrees = ntuple(i -> degree, Val(num_dims(p)))
  Quadrature(p, name, degrees; T)
end

function Quadrature(
  p::Polytope,::TensorProduct,quadratures::Vector{<:Quadrature{1}}; T::Type{<:AbstractFloat}=Float64
)
  @assert is_n_cube(p) "Tensor product quadrature rule only for n-cubes."
  D = num_dims(p)
  @assert length(quadratures) == D

  quads = map(
    q -> (map(xi -> xi[1], get_coordinates(q)), get_weights(q)),
    quadratures
  )
  coords, weights = _tensor_product(quads, T)

  names_1d = join(map(get_name, quadratures)," \n - ")
  GenericQuadrature(coords,weights,"Tensor product of 1d quadratures given by: \n "*names_1d)
end

function maxdegree(p::Polytope, ::TensorProduct)
  is_n_cube(p) ? Inf : 0
end

# Low level constructor

function _tensor_product_legendre(degrees; T::Type{<:AbstractFloat}=Float64)
  npoints = [_npoints_from_degree(degree) for degree in degrees]
  quads = [gauss(T, npoint) for npoint in npoints]
  # Geometrical map from (-1, +1) to (0, 1)
  for quad in quads
    quad[1] .+= 1
    quad[1] ./= 2
    quad[2] ./= 2
  end
  coords, weights = _tensor_product(quads, T)
  coords, weights
end

