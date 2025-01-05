
struct PolytopalQuadrature{D,T,P<:GeneralPolytope{D},Q<:Quadrature{D,T}} <: Quadrature{D,T}
  poly :: P
  quad :: Q
end

function Quadrature(p::GeneralPolytope{D},args...;kwargs...) where D
  simplex = simplex_polytope(Val(D))
  quad = Quadrature(simplex,args...;kwargs...)
  PolytopalQuadrature(p,quad)
end

# Quadrature API

# TODO: This is awful and needs to be optimised for performance. No allocations should be 
# required. We also have to optimise the mesh-wise version of these for CellQuadratures, i.e 
# 
# function Quadrature(trian::PolytopalGrid,args...;kwargs...)
#   [...]
# end

function get_coordinates(q::PolytopalQuadrature)
  cmaps = simplexify_cellmaps(q.poly)
  ref_coords = get_coordinates(q.quad)
  vcat([evaluate(cmap,ref_coords) for cmap in cmaps]...)
end

function get_weights(q::PolytopalQuadrature)
  cmaps = simplexify_cellmaps(q.poly)
  jacobians = map(gradient,cmaps)
  ref_coords = get_coordinates(q.quad)
  ref_weights = get_weights(q.quad)
  cell_measures = map(jacobians -> sum(map(meas,evaluate(jacobians,ref_coords)) .* ref_weights), jacobians)
  vcat([ref_weights .* dV for dV in cell_measures]...)
end

function get_name(q::PolytopalQuadrature{D,T,P,Q}) where {D,T,P,Q}
  "PolytopalQuadrature{$(D),$(T),$(P),$(Q)}"
end

# Utils

function simplexify_cellmaps(p::GeneralPolytope{D}) where D
  @assert !isopen(p)
  X, T = simplexify_interior(p)
  cmaps = map(T) do Tk
    Xk = X[Tk]
    origin = Xk[1]
    gradient = TensorValues.tensor_from_columns(Tuple(xk-origin for xk in Xk[2:end]))
    AffineField(gradient,origin)
  end
  cmaps
end
