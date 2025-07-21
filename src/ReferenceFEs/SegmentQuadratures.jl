
# Gauss-Jacobi and Gauss-Legendre quadratures

"""
    gauss_jacobi_quadrature(order,alpha,beta;T::Type{<:AbstractFloat}=Float64,a=zero(T),b=one(T)) -> x, w

Computes the Gauss-Jacobi quadrature rule of order `order` on the reference segment `[a,b]`
with parameters `alpha` and `beta`. Gauss-Legendre quadratures are obtained for `alpha = beta = 0`.

Uses the package `FastGaussQuadratures.jl`, which only allows Float64. If any other type 
is used, the Float64 values are converted to the requested type.
"""
function gauss_jacobi_quadrature(
  order,alpha,beta;T::Type{<:AbstractFloat}=Float64,a=zero(T),b=one(T)
)
  n = _npoints_from_order(order)
  if iszero(alpha) && iszero(beta)
    x, w = gausslegendre(n)
  else
    x, w = gaussjacobi(n,alpha,beta)
  end
  _map_segment(a,b,x,w;T)
end

"""
    gauss_legendre_quadrature(order;T::Type{<:AbstractFloat}=Float64,a=zero(T),b=one(T)) -> x, w

Computes the Gauss-Legendre quadrature rule of order `order` on the reference segment `[a,b]`.
Uses the package `QuadGK.jl`, which allows for arbitrary precision natively.
"""
function gauss_legendre_quadrature(
  order;T::Type{<:AbstractFloat}=Float64,a=zero(T),b=one(T)
)
  n = _npoints_from_order(order)
  x, w = gauss(T, n, a, b)
  x, w
end

# Utils

# Transforms a 1-D quadrature from `[-1,1]` to `[a,b]`, with `a<b`.
function _map_segment(a,b,points,weights;T::Type{<:AbstractFloat}=Float64)
  @check a < b
  points_ab = Vector{T}(undef,length(points))
  weights_ab = Vector{T}(undef,length(weights))
  points_ab .= 0.5*(b-a)*points .+ 0.5*(a+b)
  weights_ab .= 0.5*(b-a)*weights
  return points_ab, weights_ab
end

function _npoints_from_order(order)
  ceil(Int, (order + 1.0) / 2.0 )
end

# Rational quadratures

"""
    rational_gauss_legendre_quadrature(order::Integer; T::Type{<:Rational} = Rational{BigInt}, sample::Function = sample_quasi_chebyshev_3)

Computes rational Gauss-Legendre quadratures of degree `order` on the reference segment `[0,1]`.

The type of rational number if given by `T`, which defaults to `Rational{BigInt}`. The function
`sample` is used to sample the abscissae. The default is `sample_quasi_chebyshev_3`, which
is a quasi-Chebyshev sampling that keeps the weights positive for degree <= 150.
"""
function rational_gauss_legendre_quadrature(
  order::Integer; T::Type{<:Rational} = Rational{BigInt}, sample::Function = sample_quasi_chebyshev_3
)
  n = order + 1

  # Sample abscissae
  x = Vector{T}(undef, n)
  for i in 1:n
    x[i] = sample(i, n)
  end

  # Assemble Vandermonde matrix
  A = zeros(T, (n, n))
  for j in 1:n
    xj = x[j]
    p = 1 // 1
    for i in 1:n
      A[i, j] = p
      p *= xj
    end
  end

  # Assemble integrals
  b = zeros(T, n)
  for i in 1:n
    b[i] = 1 // i
  end

  # Solve for weights
  w = A \ b

  return x, w
end

# Keeps weight positivity until degree 5
sample_uniform_open(i, n) = (2 * i - 1) // (2 * n)

# Keeps weight positivity until degree 7
function sample_uniform_closed(i, n)
  if n == 1
    1 // 2
  else
    (i - 1) // (n - 1)
  end
end

# Order 2 approximation of Chebyshev nodes
# Keeps weight positivity until degree 54
function sample_quasi_chebyshev_2(i, n)
  x = (2 * i - 1) // (2 * n)

  if x < 1 / 2
    x = 2 * x^2
  else
    x = 1 - 2 * (x - 1)^2
  end

  x
end

# Order 3 approximation of Chebyshev nodes
# Keeps weight positivity until degree >= 150
function sample_quasi_chebyshev_3(i, n)
  x = (2 * i - 1) // (2 * n)

  b = 2733 // 1000 # Magic constant
  if x < 1 / 2
    x = (2 * (2 - b) * x + b) * x^2
  else
    x = 1 - (2 * (2 - b) * (1 - x) + b) * (1 - x)^2
  end

  x
end
