"""
    struct DuffyQuadrature{D,T} <: Quadrature{D,T}
      coordinates::Vector{Point{D,T}}
      weights::Vector{T}
    end

Duffy quadrature for simplices in [0,1]^D
"""
struct DuffyQuadrature{D,T} <: Quadrature{D,T}
  coordinates::Vector{Point{D,T}}
  weights::Vector{T}
end

get_coordinates(q::DuffyQuadrature) = q.coordinates

get_weights(q::DuffyQuadrature) = q.weights


"""
    DuffyQuadrature{D}(degree::Integer) where D
"""
function DuffyQuadrature{D}(order::Integer) where D
  x,w = _duffy_quad_data(order,D)
  DuffyQuadrature(x,w)
end

function _duffy_quad_data(order::Integer,D::Int)

  beta = 0
  dim_to_quad_1d = [
    _gauss_jacobi_in_0_to_1(order,(D-1)-(d-1),beta) for d in 1:(D-1) ]

  quad_1d = _gauss_legendre_in_0_to_1(order)
  push!(dim_to_quad_1d,quad_1d)

  x_pos = 1
  w_pos = 2
  dim_to_xs_1d = [quad_1d[x_pos] for quad_1d in dim_to_quad_1d]
  dim_to_ws_1d = [quad_1d[w_pos] for quad_1d in dim_to_quad_1d]

  a = 0.5
  for d in (D-1):-1:1
    ws_1d = dim_to_ws_1d[d]
    ws_1d[:] *= a
    a *= 0.5
  end

  x,w = _tensor_product_duffy(dim_to_xs_1d,dim_to_ws_1d)

  (_duffy_map.(x),w)

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

function _gauss_jacobi_in_0_to_1(order,alpha,beta)
  n = _npoints_from_order(order)
  x,w = gaussjacobi(n,alpha,beta)
  _map_to(0,1,x,w)
end

function _gauss_legendre_in_0_to_1(order)
  n = _npoints_from_order(order)
  x,w = gausslegendre(n)
  _map_to(0,1,x,w)
end

# Transforms a 1-D quadrature from `[-1,1]` to `[a,b]`, with `a<b`.
function _map_to(a,b,points,weights)
  points_ab = 0.5*(b-a)*points .+ 0.5*(a+b)
  weights_ab = 0.5*(b-a)*weights
  (points_ab, weights_ab)
end

function _npoints_from_order(order)
  ceil(Int, (order + 1.0) / 2.0 )
end

function _tensor_product_duffy(
  dim_to_xs_1d::Vector{Vector{T}},
  dim_to_ws_1d::Vector{Vector{W}}) where {T,W}

  D = length(dim_to_ws_1d)
  @assert D == length(dim_to_xs_1d)
  dim_to_n = [length(ws_1d) for ws_1d in dim_to_ws_1d]
  n = prod(dim_to_n)
  xs = zeros(Point{D,T},n)
  ws = zeros(W,n)
  cis = CartesianIndices(tuple(dim_to_n...))
  m = zero(Mutable(Point{D,T}))
  _tensor_product_duffy!(xs,ws,dim_to_xs_1d,dim_to_ws_1d,cis,m)
  (xs,ws)
end

function _tensor_product_duffy!(
  xs,ws,dim_to_xs_1d,dim_to_ws_1d,cis::CartesianIndices{D},m) where D
  k = 1
  for ci in cis
    w = 1.0
    for d in 1:D
      xs_1d = dim_to_xs_1d[d]
      ws_1d = dim_to_ws_1d[d]
      i = ci[d]
      xi = xs_1d[i]
      wi = ws_1d[i]
      w *= wi
      m[d] = xi
    end
    xs[k] = m
    ws[k] = w
    k += 1
  end
end
