
struct PolytopalQuadrature{D,T,P<:GeneralPolytope{D},Q<:Quadrature{D,T}} <: Quadrature{D,T}
  poly :: P
  quad :: Q
  conn :: Vector{Vector{Int32}}
end

function Quadrature(poly::GeneralPolytope{D},args...;kwargs...) where D
  conn, simplex = simplexify(poly)
  quad = Quadrature(simplex,args...;kwargs...)
  PolytopalQuadrature(poly,quad,conn)
end

function Quadrature(poly::GeneralPolytope{D},quad::Quadrature{D}) where D
  conn, simplex = simplexify(poly)
  @check get_polytope(quad) == simplex
  PolytopalQuadrature(poly,quad,conn)
end

# Quadrature API

function get_name(q::PolytopalQuadrature{D,T,P,Q}) where {D,T,P,Q}
  "PolytopalQuadrature{$(D),$(T),$(P),$(Q)}"
end

get_coordinates(q::PolytopalQuadrature) = evaluate(PQuadCoordsMap(),q)

function Arrays.lazy_map(::typeof(get_coordinates),quads::AbstractVector{<:PolytopalQuadrature})
  lazy_map(PQuadCoordsMap(),quads)
end

get_weights(q::PolytopalQuadrature) = evaluate(PQuadWeightsMap(),q)

function Arrays.lazy_map(::typeof(get_weights),quads::AbstractVector{<:PolytopalQuadrature})
  lazy_map(PQuadWeightsMap(),quads)
end

# Private functions

struct PQuadCoordsMap <: Map end

function Arrays.return_cache(::PQuadCoordsMap, q::PolytopalQuadrature)
  conn = q.conn
  coords = get_vertex_coordinates(q.poly)
  cmap = affine_map(Tuple(coords[first(conn)]))
  
  x = get_coordinates(q.quad)
  cmap_cache = return_cache(cmap,x)
  y = evaluate!(cmap_cache,cmap,x)
  y_cache = CachedArray(similar(y,(length(x)*length(conn))))

  return y_cache, cmap_cache
end

function Arrays.evaluate!(cache,::PQuadCoordsMap, q::PolytopalQuadrature)
  y_cache, cmap_cache = cache
  coords = get_vertex_coordinates(q.poly)
  conn = q.conn

  x = get_coordinates(q.quad)
  setsize!(y_cache,(length(x)*length(conn),))
  y = y_cache.array

  nx = length(x)
  for (k,verts) in enumerate(conn)
    cmap = affine_map(Tuple(coords[verts]))
    y[(k-1)*nx+1:k*nx] .= evaluate!(cmap_cache, cmap, x)
  end
  return y
end

struct PQuadWeightsMap <: Map end

function Arrays.return_cache(::PQuadWeightsMap, q::PolytopalQuadrature)
  conn = q.conn
  w = get_weights(q.quad)
  return CachedArray(similar(w,(length(w)*length(conn))))
end

function Arrays.evaluate!(cache,::PQuadWeightsMap, q::PolytopalQuadrature)
  conn = q.conn
  coords = get_vertex_coordinates(q.poly)

  w = get_weights(q.quad)
  setsize!(cache,(length(w)*length(conn),))
  z = cache.array

  nw = length(w)
  for (k,verts) in enumerate(conn)
    dV = meas(affine_map(Tuple(coords[verts])).gradient)
    z[(k-1)*nw+1:k*nw] .= w .* dV
  end
  return z
end
