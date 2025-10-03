
struct Duffy <: QuadratureName end
const duffy = Duffy()

function Quadrature(
  p::Polytope, ::Duffy, degree::Integer; T::Type{<:AbstractFloat}=Float64
)
  msg = """\n
  `duffy` quadrature rule only for simplices.
  """
  @assert is_simplex(p) msg

  D = num_dims(p)
  coords, weights = _duffy_quad(degree, D; T)
  name = "Simplex quadrature of degree $degree obtained with Duffy transformation and tensor product of 1d Gauss-Jacobi and Gauss-Legendre rules."
  GenericQuadrature(coords, weights, name)
end

function maxdegree(p::Polytope, ::Duffy)
  is_simplex(p) ? Inf : 0
end

function _duffy_quad(degree::Integer, D::Int; T::Type{<:AbstractFloat}=Float64)
  beta = 0
  quads = [
    _gauss_jacobi_in_0_to_1(degree, (D - 1) - (d - 1), beta; T)
    for d in 1:(D-1)
  ]
  quad_1d = _gauss_legendre_in_0_to_1(degree; T)
  push!(quads, quad_1d)

  a = one(T) / 2
  for d in (D-1):-1:1
    rmul!(quads[d][2], a)
    a /= 2
  end

  _coords, weights = _tensor_product(quads, T)
  coords = _duffy_map.(_coords)
  coords, weights
end

# Duffy map from the n-cube in [0,1]^d to the n-simplex in [0,1]^d
function _duffy_map(q::Point{D,T}) where {D,T}
  m = zero(Mutable(Point{D,T}))
  m[1] = q[1]
  a = one(T)
  for i in 2:D
    a *= (1 - q[i-1])
    m[i] = a * q[i]
  end
  Point{D,T}(m)
end

_duffy_map(q::Point{1,T}) where T = q

function _gauss_jacobi_in_0_to_1(degree, alpha, beta; T::Type{<:AbstractFloat}=Float64)
  n = _npoints_from_degree(degree)
  x, w = gaussjacobi(n, alpha, beta)
  _map_to(0, 1, x, w; T)
end

function _gauss_legendre_in_0_to_1(degree; T::Type{<:AbstractFloat}=Float64)
  n = _npoints_from_degree(degree)
  x, w = gausslegendre(n)
  _map_to(0, 1, x, w; T)
end

# Transforms a 1-D quadrature from `[-1,1]` to `[a,b]`, with `a<b`.
function _map_to(a, b, points, weights; T::Type{<:AbstractFloat}=Float64)
  points_ab = Vector{T}(undef, length(points))
  weights_ab = Vector{T}(undef, length(weights))
  c = (a + b) / 2
  h = (b - a) / 2
  points_ab .= h * points .+ c
  weights_ab .= h * weights
  points_ab, weights_ab
end
