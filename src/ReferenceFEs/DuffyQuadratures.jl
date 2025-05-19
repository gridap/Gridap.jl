
"""
    struct Duffy <: QuadratureName

Duffy quadrature rule for simplices, obtained as the mapped
tensor product of 1d Gauss-Jacobi and Gauss-Legendre quadratures.

# Constructor: 

    Quadrature(p::Polytope,duffy,degree::Integer;T::Type{<:AbstractFloat}=Float64)
"""
struct Duffy <: QuadratureName end

const duffy = Duffy()

function Quadrature(p::Polytope,::Duffy,degree::Integer;T::Type{<:AbstractFloat}=Float64)
  @assert is_simplex(p)
  D = num_dims(p)
  x, w = _duffy(degree,D,T=T)
  name = "Simplex quadrature of degree $degree obtained with Duffy transformation and tensor product of 1d Gauss-Jacobi and Gauss-Legendre rules."
  GenericQuadrature(x,w,name)
end

function _duffy(degree::Integer,D::Int;T::Type{<:AbstractFloat}=Float64)
  quad_1d = [ gauss_jacobi_quadrature(degree,(D-1)-(d-1),0;T) for d in 1:D]
  coords_1d, weights_1d = map(first, quad_1d), map(last, quad_1d)
  
  a = T(0.5)
  for d in (D-1):-1:1
    weights_1d[d] .*= a
    a *= T(0.5)
  end

  coords, weights = _tensor_product(coords_1d, weights_1d)
  return _duffy_map.(coords), weights
end

# Duffy map from the n-cube in [0,1]^d to the n-simplex in [0,1]^d
function _duffy_map(q::Point{D,T}) where {D,T}
  m = zero(Mutable(Point{D,T}))
  m[1] = q[1]
  a = one(T)
  for i in 2:D
    a *= (1-q[i-1])
    m[i] = a*q[i]
  end
  Point{D,T}(m)
end

_duffy_map(q::Point{1,T}) where T = q
