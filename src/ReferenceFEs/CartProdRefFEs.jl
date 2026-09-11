# The Cartesian product V₁ × … × V_K of reference elements; see `CartProdRefFE`.
#
# Two conventions, aligned: the copy index is **last** inside a value, and factor
# `c` occupies **block `c`** of the numbering. Both follow from Gridap writing
# derivative indices first, ∇A[k,i,j] = ∂_k A_ij. Appending then commutes with ∇,
# so one rule, `outer(v, e)`, serves a value, a gradient and a Hessian alike;
# prepending would need the copy index stepped over the derivative indices, which
# no contraction of existing operations produces, and would transpose every
# gradient invisibly whenever K == D.
#
# Together they are what Gridap's own operators expect. `divergence` is `tr(∇·)`
# and `tr` traces the first two indices, so `div(A)_j = ∂ᵢ A_ij` is the divergence
# of the *columns*, i.e. the vector of the factors' divergences. And
# `_generate_dof_layout_node_major` blocks by component, so a stacked scalar
# Lagrangian element reproduces `lagrangian(VectorValue{K,T})` DoF for DoF.
#
# The price: the weak-symmetry stress space is usually written with one H(div)
# *row* per component of the divergence; here it is one per column.

# The injection ι_c(v) = v ⊗ e_c and its inverse π_c(v) = v ⋅ e_c, the two
# constant linear maps the whole construction rests on. Shared by the basis
# scatter and the moment stacking -- the same operation, which is why the DoFs
# need no type of their own. Only rank ≤ 2 factors stack: a 2-tensor-valued one
# would need a rank-4 stacked gradient, which `outer` does not build.
@inline _cp_insert(e, v) = outer(v, e)
@inline _cp_extract(e, v) = v ⋅ e

# `E` throughout is `representatives_of_componentbasis_dual(VectorValue{K})`, the
# stacking axis. Every use needs all K of it at once, so it is built whole and
# kept in whatever cache is already at hand. That basis is orthonormal, so the
# same `e` both injects and reads back.

############################################################################################
# The stacked basis

# `K` stacked bases, factor `c` in block `c`. `bases` is either one basis, of
# which the stack is `K` copies, or a `K`-tuple of them. The two cases differ in
# one place only -- a power evaluates once and scatters, a product evaluates each
# factor -- so they are one struct and two aliases, not two types.
#
# It must be a type rather than a `lazy_map` over the factors' elements: a basis
# built on a `PolynomialBasis` does not support element access -- `getindex`
# there returns the dummy `PT()` rather than a field -- so going through it would
# produce garbage silently. Only array-level evaluation is safe, so that is all
# this offers.
struct CartProdBasis{K,B} <: AbstractVector{Field}
  bases::B
end

CartProdBasis{K}(b::B) where {K,B} = CartProdBasis{K,B}(b)

const CartProdPowerBasis{K} = CartProdBasis{K,<:AbstractVector{<:Field}} # V^K
const CartProdTupleBasis{K} = CartProdBasis{K,<:Tuple} # V_1 × … × V_K

Base.IndexStyle(::Type{<:CartProdBasis}) = IndexLinear()

Base.size(b::CartProdPowerBasis{K}) where K = (K * length(b.bases),)
Base.size(b::CartProdTupleBasis) = (sum(length, b.bases),)

get_order(b::CartProdPowerBasis) = get_order(b.bases)
get_order(b::CartProdTupleBasis) = maximum(get_order, b.bases)

Base.getindex(b::CartProdPowerBasis, ::Integer) = first(b.bases)
Base.getindex(b::CartProdTupleBasis, ::Integer) = first(first(b.bases))

# values and derivatives share one scatter, differing only in what they evaluate
_cp_src(f, ::Val{0}) = f
_cp_src(f, ::Val{1}) = Broadcasting(∇)(f)
_cp_src(f, ::Val{2}) = Broadcasting(∇∇)(f)

# The cache holds the sources and their caches -- one shared for a power, one per
# factor for a product -- and `_cp_vals` turns them into the `K` value arrays
# without evaluating a power `K` times.
function _cp_prepare(b::CartProdPowerBasis, N, x)
  src = _cp_src(b.bases, N)
  return src, return_cache(src, x)
end

function _cp_prepare(b::CartProdTupleBasis, N, x)
  srcs = map(f -> _cp_src(f, N), b.bases)
  return srcs, map(f -> return_cache(f, x), srcs)
end

function _cp_vals(::CartProdPowerBasis{K}, cs, src, x) where K
  v = evaluate!(cs, src, x)
  return ntuple(i -> v, K)
end

_cp_vals(::CartProdTupleBasis, cs, srcs, x) =
  map((c, f) -> evaluate!(c, f, x), cs, srcs)

function _cp_cache(b::CartProdBasis{K}, x, N) where K
  srcs, cs = _cp_prepare(b, N, x)
  vs = _cp_vals(b, cs, srcs, x)
  E = representatives_of_componentbasis_dual(VectorValue{K,Float64})
  # every factor shares a value type, which is what makes the stack a MultiValue
  T = typeof(_cp_insert(first(E), testitem(first(vs))))
  r = CachedArray(zeros(T, size(first(vs), 1), sum(v -> size(v, 2), vs)))
  return r, cs, srcs, E
end

function _cp_eval!(cache, b::CartProdBasis{K}, x, N) where K
  r, cs, srcs, E = cache
  vs = _cp_vals(b, cs, srcs, x)
  np = size(first(vs), 1)
  setsize!(r, (np, sum(v -> size(v, 2), vs)))
  a = r.array
  off = 0
  @inbounds for c in 1:K
    v, e = vs[c], E[c]
    for j in axes(v, 2), i in 1:np
      a[i, off+j] = _cp_insert(e, v[i, j])
    end
    off += size(v, 2)
  end
  return a
end

return_cache(b::CartProdBasis, x::AbstractVector{<:Point}) = _cp_cache(b, x, Val(0))

evaluate!(cache, b::CartProdBasis, x::AbstractVector{<:Point}) =
  _cp_eval!(cache, b, x, Val(0))

return_cache(
  fg::FieldGradientArray{N,<:CartProdBasis}, x::AbstractVector{<:Point}
) where N = _cp_cache(fg.fa, x, Val(N))

evaluate!(
  cache, fg::FieldGradientArray{N,<:CartProdBasis}, x::AbstractVector{<:Point}
) where N = _cp_eval!(cache, fg.fa, x, Val(N))

############################################################################################
# The stacked DoF basis

# The primitive is one slot, not the whole stack: the DoFs of `b` acting on slice
# `c` of a `K`-fold stacked value,
#
#   σ_{i,c}(u) = σ_i(π_c u),    π_c u = u ⋅ e_c,
#
# which is an ordinary DoF basis of the same kind as `b` -- the slicing is
# absorbed into the moment weights, possible because a moment DoF is linear in
# the field. Nothing here has to differentiate, cache, or be evaluated
# field-by-field, and the DoF indices within a slot are the atom's own, so face
# ownership carries over untouched.
#
# The stack is then the `vcat` of the `K` slots, which is what puts copy `c` in
# block `c`. It has to be a `vcat`: a single `MomentBasedDofBasis` emits its DoFs
# *face-major* (see its `evaluate!`), so one stacked basis could only interleave
# the copies within each face, never block them across the element.
function _cp_slot_dofs(b::MomentBasedDofBasis, c::Integer, ::Val{K}) where K
  e = representatives_of_componentbasis_dual(VectorValue{K,Float64})[c]
  fm = get_face_moments(b)
  W = typeof(_cp_insert(e, zero(eltype(eltype(fm)))))
  moments = map(M -> W[_cp_insert(e, m) for m in M], fm)
  return MomentBasedDofBasis(
    get_nodes(b), moments, get_face_nodes_dofs(b), get_face_own_moments(b), b.operator
  )
end

_cp_slot_dofs(b::ConcatenatedDofVector, c::Integer, K::Val) =
  vcat((_cp_slot_dofs(a, c, K) for a in b.args)...)

_cp_slot_dofs(b::LinearCombinationDofVector, c::Integer, K::Val) =
  linear_combination(b.values, _cp_slot_dofs(b.predofs, c, K))

# A nodal basis has no slot form: `LagrangianDofBasis` assigns a DoF to *every*
# component of its value type, so it cannot leave the other `K-1` slots empty.
# Its whole stack is built at once below, and one slot is a slice of that --
# contiguous, because the copies are blocked.
function _cp_slot_dofs(b::LagrangianDofBasis, c::Integer, K::Val)
  n = length(b.dof_to_node)
  restrict(_cp_stack_dofs(b, K), ((c-1)*n+1):(c*n))
end

_cp_slot_dofs(b, c, ::Val) = @notimplemented """\n
Do not know how to slice a $(typeof(b)); add a `_cp_slot_dofs` method for it.
"""

# `K` different DoF bases, factor `i` on slot `i`.
_cp_stack_dofs(bs::Tuple, K::Val{NK}) where NK =
  vcat(ntuple(c -> _cp_slot_dofs(bs[c], c, K), NK)...)

# `K` copies of one, which is the same thing with every factor equal
_cp_stack_dofs(b, K::Val{NK}) where NK = _cp_stack_dofs(ntuple(i -> b, NK), K)

# ... except for a nodal basis, which stacks whole with no slicing at all. For an
# atom in Gridap's canonical layout the result is exactly
# `LagrangianDofBasis(W, b.nodes)`; it is built explicitly so that a
# non-canonical atom is stacked correctly too.
#
# With `d = num_indep_components(V)`, DoF `a` of the atom reads component
# `j = dof_to_comp[a]` at node `dof_to_node[a]`, and copy `c` of it reads
# component `j + d*(c-1)` of the stacked value -- the linear independent-component
# index of `(j,c)`, since the copy index goes last and `MultiValue`s are column
# major.
function _cp_stack_dofs(b::LagrangianDofBasis{P,V}, ::Val{K}) where {P,V,K}
  ndofs = length(b.dof_to_node)
  nnodes = length(b.nodes)
  d = num_indep_components(V)

  e = zero(VectorValue{K,Float64})
  W = change_eltype(typeof(_cp_insert(e, zero(V))), Int)
  ncomps = num_indep_components(W)

  dof_to_node = zeros(Int, K * ndofs)
  dof_to_comp = zeros(Int, K * ndofs)
  m = zeros(Int, nnodes, ncomps)
  for a in 1:ndofs, c in 1:K
    dof = (c - 1) * ndofs + a
    comp = b.dof_to_comp[a] + d * (c - 1)
    dof_to_node[dof] = b.dof_to_node[a]
    dof_to_comp[dof] = comp
    m[b.dof_to_node[a], comp] = dof
  end
  node_and_comp_to_dof = [W(ntuple(l -> m[node, l], ncomps)) for node in 1:nnodes]

  LagrangianDofBasis(b.nodes, dof_to_node, dof_to_comp, node_and_comp_to_dof)
end

############################################################################################
# The push-forward

# Factor `i`'s push-forward applied to slice `i`,
#
#   φ = Σ_i ι_i( PF_i(π_i φ̂) ),
#
# with `PFS` the `Tuple` type of the factors' maps. For `V^K` they are `K` copies
# of one.
#
# This is *not* one of the double Piola maps. Stacking K contravariant Piola maps
# gives φ = det(J)⁻¹ J φ̂, contravariant on the first index only, whereas
# `DoubleContraVariantPiolaMap` is det(J)⁻² J φ̂ Jᵀ. So Hellan-Herrmann-Johnson is
# not obtainable by wrapping Raviart-Thomas, and should not be: the column-wise
# tensor is a different element with a different conformity.
struct CartProdPushforward{K,PFS<:Tuple} <: Pushforward end

# `IdentityPiolaMap` is short-circuited rather than evaluated by Gridap, so a
# stack that mixes it with a Piola map cannot take that branch and must apply it
# explicitly. It is the identity.
@inline _cp_apply(::IdentityPiolaMap, v, Jt) = v
@inline _cp_apply(::InversePushforward{IdentityPiolaMap}, v, Jt) = v
@inline _cp_apply(pf, v, Jt) = evaluate(pf, v, Jt)

_cp_maps(::Type{PFS}) where PFS<:Tuple =
  ntuple(i -> fieldtype(PFS, i)(), fieldcount(PFS))

_cp_blockwise(pfs, v, Jt, E) =
  sum(_cp_insert(E[i], _cp_apply(pfs[i], _cp_extract(E[i], v), Jt)) for i in eachindex(E))

# `OperationField` builds the cache of its operation once and hands it back to
# every `evaluate!`, so these are built once per cell, not per value.
return_cache(::CartProdPushforward{K,PFS}, ::Number, ::Number) where {K,PFS} =
  (representatives_of_componentbasis_dual(VectorValue{K,Float64}), _cp_maps(PFS))

evaluate!(
  cache, ::CartProdPushforward{K,PFS}, v_ref::Number, Jt::Number
) where {K,PFS} = _cp_blockwise(cache[2], v_ref, Jt, cache[1])

return_cache(
  ::InversePushforward{CartProdPushforward{K,PFS}}, ::Number, ::Number
) where {K,PFS} = (representatives_of_componentbasis_dual(VectorValue{K,Float64}),
                   map(inverse_map, _cp_maps(PFS)))

evaluate!(
  cache, ::InversePushforward{CartProdPushforward{K,PFS}}, v_phys::Number, Jt::Number
) where {K,PFS} = _cp_blockwise(cache[2], v_phys, Jt, cache[1])

############################################################################################
# The reference FE

"""
    struct CartProd{K,NS} <: ReferenceFEName end

Reference FE name for `K` stacked elements. `NS` is the factors' own name type
when they are all the same element (`V^K`), and the `Tuple` of their name types
otherwise (`V₁ × … × V_K`). See [`CartProdRefFE`](@ref).
"""
struct CartProd{K,NS} <: ReferenceFEName end

_cp_name_maps(::Val{K}, ::Type{N}) where {K,N<:ReferenceFEName} =
  ntuple(i -> Pushforward(N), K)

_cp_name_maps(::Val{K}, ::Type{NS}) where {K,NS<:Tuple} =
  ntuple(i -> Pushforward(fieldtype(NS, i)), K)

function Pushforward(::Type{CartProd{K,NS}}) where {K,NS}
  pfs = _cp_name_maps(Val(K), NS)
  # K copies of the identity is the identity, and Gridap short-circuits that map
  # rather than evaluating it, so it must be returned unwrapped.
  all(pf -> isa(pf, IdentityPiolaMap), pfs) && return IdentityPiolaMap()
  return CartProdPushforward{K,Tuple{map(typeof, pfs)...}}()
end

"""
    CartProdRefFE(reffe::ReferenceFE, K::Integer)
    CartProdRefFE(reffes::ReferenceFE...)

`K` stacked elements, along a new last index of the value type — the reference FE
of the Cartesian product space `V₁ × … × V_K`, or of the power `V^K` when a
single element and a count are given. Values are

    Float64        ⟹  VectorValue{K}
    VectorValue{d} ⟹  TensorValue{d,K}

so that factor `c` of a vector-valued atom is the `c`-th *column*. Scalar factors
give vector-valued Morley, Argyris or Hermite, and a displacement whose
components have different regularity; vector-valued ones give the column-wise
H(div) tensors of elasticity with weakly imposed symmetry [Arnold, Falk &
Winther, Numer. Math. 92 (2002) 401].

All factors must share a polytope, a `Conformity` and a value type; the last is
what makes the result a single `MultiValue` rather than a `MultiFieldFESpace`.

The DoFs, shape functions, prebasis and face ownership are blocked by factor,
with item `i` of factor `c` at index `offset(c) + i`. The conformity and polytope
are inherited, and the factors are kept as the metadata, which is how
`compute_cell_bases_changes` recovers them.

## Examples

    CartProdRefFE(MorleyRefFE(Float64, TRI), 2)        # vector-valued Morley
    CartProdRefFE(ReferenceFE(TRI, raviart_thomas, Float64, 1), 2)
                                                       # column-wise H(div) tensor
    CartProdRefFE(ReferenceFE(TRI, lagrangian, Float64, 2),
                  ReferenceFE(TRI, lagrangian, Float64, 2),
                  ArgyrisRefFE(Float64, TRI))          # Kirchhoff--Love shell
"""
function CartProdRefFE(reffes::ReferenceFE...)
  K = length(reffes)
  @notimplementedif K < 1 "Need at least one factor."
  @notimplementedif any(r -> get_polytope(r) != get_polytope(first(reffes)), reffes) """\n
  All factors must share a polytope.
  """
  @notimplementedif any(r -> Conformity(r) != Conformity(first(reffes)), reffes) """\n
  All factors must share a conformity, got $(map(Conformity, reffes)).
  """
  @notimplementedif any(r -> _cp_value_type(r) != _cp_value_type(first(reffes)), reffes) """\n
  All factors must share a value type, got a mix. Use a `MultiFieldFESpace`.
  """
  NS = Tuple{map(r -> typeof(get_name(r)), reffes)...}
  _cp_reffe(CartProd{K,NS}, reffes,
            CartProdBasis{K}(map(get_prebasis, reffes)),
            CartProdBasis{K}(map(get_shapefuns, reffes)))
end

function CartProdRefFE(reffe::ReferenceFE, ::Val{K}) where K
  @notimplementedif K < 1 "Need at least one copy, got $K."
  N = typeof(get_name(reffe))
  _cp_reffe(CartProd{K,N}, ntuple(i -> reffe, K),
            CartProdBasis{K}(get_prebasis(reffe)),
            CartProdBasis{K}(get_shapefuns(reffe)))
end

CartProdRefFE(reffe::ReferenceFE, K::Integer) = CartProdRefFE(reffe, Val(K))

_cp_value_type(reffe) = typeof(testitem(
  evaluate(get_shapefuns(reffe), get_vertex_coordinates(get_polytope(reffe)))
))

# The factors never mix, so everything is blocked: the DoFs are the `vcat` of the
# slots, and a face owns each factor's own DoFs shifted past the factors before
# it.
function _cp_reffe(::Type{Name}, reffes, prebasis, shapefuns) where Name
  K = length(reffes)
  p = get_polytope(first(reffes))
  conf = Conformity(first(reffes))
  offs = _cp_offsets(map(num_dofs, reffes))

  dofs = _cp_stack_dofs(map(get_dof_basis, reffes), Val(K))
  owns = map(r -> get_face_own_dofs(r, conf), reffes)
  face_own_dofs = [
    Int[offs[c] + d for c in 1:K for d in owns[c][f]] for f in eachindex(first(owns))
  ]

  return GenericRefFE{Name}(
    sum(num_dofs, reffes), p, prebasis, dofs, conf, reffes, face_own_dofs, shapefuns
  )
end

# `_cp_offsets(ns)[c]` is the number of items before factor `c`
_cp_offsets(ns) = (0, cumsum(ns)[1:end-1]...)

# The factors' permutations, blocked. The offsets are per face: factor `c`'s own
# DoFs on a face sit after those of the factors before it *on that same face*.
# `INVALID_PERM` propagates.
function get_face_own_dofs_permutations(
  reffe::GenericRefFE{CartProd{K,NS}}, conf::Conformity
) where {K,NS}
  reffes = get_metadata(reffe)
  perms = map(r -> get_face_own_dofs_permutations(r, conf), reffes)
  owns = map(r -> get_face_own_dofs(r, conf), reffes)
  map(eachindex(first(perms))) do f
    offs = _cp_offsets(map(o -> length(o[f]), owns))
    map(eachindex(first(perms)[f])) do pindex
      Int[q == INVALID_PERM ? INVALID_PERM : offs[c] + q
          for c in 1:K for q in perms[c][f][pindex]]
    end
  end
end
