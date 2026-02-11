
struct TensorProduct <: QuadratureName end
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
