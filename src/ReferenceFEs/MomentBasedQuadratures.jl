# Moment-based quadratures.
#
# We assume that:
#  - Both the quadrature and the basis within the map are defined on the reference space
#    of the polytope. Only QUAD and HEX are supported.
#  - The inputs are defined on the physical space, since the quadrature is
#    suitable to integrate polynomials on the polytope.
# Kronecker is a flag to indicate whether the basis has a Kronecker property
# on the nodes of the quadrature. If true, we can skip the inversion of the lhs matrix.

"""
    struct MomentBasedQuadratureMap{P,Q,B} <: Map

Map that, given a polytope `poly` and a quadrature `quad` able to integrate polynomials on it,
returns a quadrature with the nodes of a Lagrangian basis on the bounding box of `poly`
and weights matching the moments of that basis computed with `quad`. `poly` can also be
a vector of points (e.g. the physical coordinates of a cell), in which case their bounding box is used.

    evaluate(k::MomentBasedQuadratureMap, poly, quad::Quadrature) -> GenericQuadrature
"""
struct MomentBasedQuadratureMap{P,Q,B} <: Map
  poly  :: P
  nodes :: Q
  basis :: B
  kronecker :: Bool
end

"""
    MomentBasedQuadratureMap(D::Integer, degree; T=Float64, space=:Q)

Returns a [`MomentBasedQuadratureMap`](@ref) in dimension `D`, exact for the polynomial
space `space ∈ (:P,:Q)` of order `degree`.
"""
function MomentBasedQuadratureMap(
  D::Integer, degree;
  T::Type = Float64, space::Symbol = :Q
)
  if space == :P
    ref_poly = ifelse(D==2,TRI,TET)
    filter = Polynomials._p_filter
  elseif space == :Q
    ref_poly = ifelse(D==2,QUAD,HEX)
    filter = Polynomials._q_filter
  else
    error("Space $space not supported.")
  end
  ref_dofs = LagrangianDofBasis(T,ref_poly,degree)
  ref_prebasis = Polynomials.JacobiPolynomialBasis{D}(T,degree,filter)
  ref_basis = compute_shapefuns(ref_dofs,ref_prebasis)
  return MomentBasedQuadratureMap(ref_poly, ref_dofs.nodes, ref_basis, true)
end

function Arrays.return_cache(k::MomentBasedQuadratureMap, poly, quad::Quadrature)
  cmap = get_bounding_box_cmap(poly)
  invcmap = inverse_map(cmap)

  cmap_cache = return_cache(cmap, k.nodes)

  xq_phys, wq = get_coordinates(quad), get_weights(quad)
  invcmap_cache = return_cache(invcmap, xq_phys)

  xq_ref = evaluate!(invcmap_cache, invcmap, xq_phys)
  i_cache = return_cache(integrate, k.basis, xq_ref, wq)

  M = (!k.kronecker) ? lu(evaluate(k.basis, k.nodes)) : nothing

  return cmap_cache, invcmap_cache, i_cache, M
end

function Arrays.evaluate!(cache, k::MomentBasedQuadratureMap, poly, quad::Quadrature)
  cmap_cache, invcmap_cache, i_cache, M = cache

  cmap = get_bounding_box_cmap(poly)
  invcmap = inverse_map(cmap)

  x_phys = evaluate!(cmap_cache, cmap, k.nodes)

  xq_phys, wq = get_coordinates(quad), get_weights(quad)
  xq_ref = evaluate!(invcmap_cache, invcmap, xq_phys)
  w_phys = evaluate!(i_cache, integrate, k.basis, xq_ref, wq)

  if !k.kronecker
    ldiv!(M, w_phys)
  end

  name = "Moment-Based quadrature"
  return GenericQuadrature(x_phys, w_phys, name)
end

function get_bounding_box_cmap(poly)
  pmin, pmax = get_bounding_box(poly)
  O = pmin
  G = diagonal_tensor(pmax - pmin)
  return affine_map(G, O)
end

"""
    MomentBasedQuadrature(p::Polytope{D}, degree; T=Float64, space=:Q)

Returns a moment-based quadrature on `p`, exact for the polynomial space `space ∈ (:P,:Q)`
of order `degree`. The moments are computed with `Quadrature(p,degree;T=T)`.
"""
function MomentBasedQuadrature(
  p::Polytope{Dc}, degree;
  T::Type = Float64, space::Symbol = :Q
) where Dc
  k = MomentBasedQuadratureMap(Dc, degree; T=T, space=space)
  quad = Quadrature(p, degree; T=T)
  return evaluate(k,p,quad)
end
