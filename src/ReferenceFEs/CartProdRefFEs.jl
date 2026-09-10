# K independent copies of a reference FE, stacked into a bigger value type.
#
# Given a reference FE whose shape functions take values in `S`, this builds the
# one whose shape functions take values in `T`, the K-fold stack of `S` along a
# new *last* index:
#
#   S = Float64            ⟹  T = VectorValue{K}
#   S = VectorValue{d}     ⟹  T = TensorValue{d,K}
#
# The FE space is the K-fold Cartesian product V ⊕ … ⊕ V, which is why the name
# matches `Polynomials.CartProdPolyBasis` -- that type is this same operation one
# module down. The atom `S` may be scalar, giving vector-valued
# Morley/Argyris/Hermite for shells and strain-gradient elasticity; or it may be
# vector-valued, giving the column-wise H(div) tensors of elasticity with weakly
# imposed symmetry [Arnold, Falk & Winther]. Nothing here cares which, only that
# T decomposes as S^K linearly.
#
# WHAT MAKES IT WORK, AND WHY THERE IS ONLY ONE NEW TYPE
#
# Everything rests on one constant linear map, the injection ι_c(s) = s ⊗ e_c,
# and on the fact that it commutes with differentiation. Two consequences: the
# generalized Vandermonde is `kron(V, I_K)`, since σ_{i,c}(Φ_{j,c'}) = δ_{cc'}
# σ_i(φ_j) -- the copies never see each other -- and so is the change of basis,
# provided the push-forward acts blockwise, which is what `CartProdPushforward`
# is for.
#
# The DoFs need *no* type at all. A moment DoF is σ(u) = Σᵢ mᵢ ⊙ u(xᵢ), so
# slicing the field is the same as stacking the weight,
#
#   σᵢ(π_c u) = Σ (mᵢ ⊗ e_c) ⊙ u(xᵢ),
#
# and the stacked DoF basis is an ordinary `MomentBasedDofBasis`. The other DoF
# kinds recurse: a `ConcatenatedDofVector` is the `vcat` of its stacked blocks,
# and a `LinearCombinationDofVector` carries `kron` on its coefficient matrix. So
# the only thing that must be wrapped is the *basis*, and only because a basis
# has to be evaluated array-wise.
#
# WHERE THE COPY INDEX GOES: LAST, ALWAYS
#
# A moment must live in the same space as the value it pairs with, and Gridap
# writes derivative indices *first*: ∇A[k,i,j] = ∂_k A_ij. An injection that
# appends the copy index therefore commutes with differentiation -- appending on
# the right never collides with prepending on the left -- and one rule covers
# every atom and every operator:
#
#   atom    operator   value layout          insertion
#   scalar  --         VectorValue{K}        outer(v, e)
#   scalar  ∇          TensorValue{D,K}      outer(v, e)
#   scalar  ∇∇         ThirdOrder{D,D,K}     outer(v, e)
#   vector  --         TensorValue{d,K}      outer(v, e)
#   vector  ∇          ThirdOrder{D,d,K}     outer(v, e)
#
# Prepending instead would need the copy index stepped over the derivative
# indices, a middle insertion no contraction of existing operations produces. It
# would also transpose every gradient, invisibly whenever K == D.
#
# This is the convention `Polynomials.CartProdPolyBasis` already uses:
# `_cartprod_set_derivative!` builds its gradient components as `∇i ⊗ vj`.
#
# It is also the one Gridap's operators expect. `divergence(f) = tr(∇f)` and `tr`
# of a third-order tensor traces its *first two* indices, so Gridap's tensor
# divergence is div(A)_j = ∂ᵢ A_ij, the divergence of the columns of A. Here the
# copies are the columns and the stacked gradient is indexed (k,j,c), so `tr`
# contracts the derivative index against the atom's own and returns one
# divergence per copy: `divergence(σₕ)` is exactly the vector of the copies'
# divergences. The cost is that the weak-symmetry stress space is conventionally
# written with one H(div) row per component of the divergence, and here it is one
# per column -- weak forms are transposed accordingly.
#
# Only rank ≤ 2 atoms stack: an atom valued in a 2-tensor would need a rank-4
# stacked gradient, which `outer` does not build and says so itself.
#
# ORDERING CONVENTION
#
# Copy `c` of base DoF/basis function `i` sits at index `K*(i-1) + c`, i.e. the
# copy index runs *fastest*. This matches `_cartprod_set_value!` in
# `Polynomials`, which scatters a term's components consecutively, so the DoF
# ordering and the polynomial-basis ordering agree. Getting these two out of step
# transposes `kron(P, I_K)` into `kron(I_K, P)`: everything is subtly wrong while
# the dimensions still check out. Note this is independent of the layout question
# above -- one is the order in which basis functions and DoFs are numbered, the
# other the arrangement of indices inside a single value.

# The injection ι_c(v) = v ⊗ e_c and its inverse π_c(v) = v ⋅ e_c, the two
# constant linear maps the whole construction rests on. No case analysis is
# needed: appending the copy index commutes with Gridap's prepending of
# derivative indices, so one rule serves a value, a gradient and a Hessian alike.
# Shared by the basis scatter and the moment stacking -- the same operation,
# which is why the DoFs need no type of their own.
@inline _cp_insert(e, v) = outer(v, e)
@inline _cp_extract(e, v) = v ⋅ e

# `E` throughout is `representatives_of_componentbasis_dual(VectorValue{K})`, the
# stacking axis. Every use needs all K of it at once, so it is built whole and
# kept in whatever cache is already at hand. That basis is orthonormal, so the
# same `e` both injects and reads back.

############################################################################################
# The stacked basis -- the only new type

# `K` copies of every field of `base`, with copy `c` of base function `j` at
# index `K*(j-1) + c`.
#
# It must be a type rather than a `lazy_map` over the elements of `base`: a basis
# built on a `PolynomialBasis` does not support element access -- `getindex`
# there returns the dummy `PT()` rather than a field -- so going through it would
# produce garbage silently. Only array-level evaluation is safe, so that is all
# this offers.
struct CartProdBasis{K,B} <: AbstractVector{Field}
  base::B
end

CartProdBasis{K}(b::B) where {K,B} = CartProdBasis{K,B}(b)

Base.size(b::CartProdBasis{K}) where K = (K * length(b.base),)
Base.IndexStyle(::Type{<:CartProdBasis}) = IndexLinear()

# A basis is atomic here: nothing should ever *evaluate* an element. Gridap still
# needs one to name types with -- `FieldGradientArray{N}(f)` asks for it -- so
# these return a base field as a type witness, exactly as `PolynomialBasis`
# returns `PT()`. The value is nonsense; never evaluate it.
Gridap.Arrays.testitem(b::CartProdBasis) = testitem(b.base)
Base.getindex(b::CartProdBasis, ::Integer) = testitem(b)

get_order(b::CartProdBasis) = get_order(b.base)

# values and derivatives share one scatter, differing only in what they evaluate
_cp_src(b::CartProdBasis, ::Val{0}) = b.base
_cp_src(b::CartProdBasis, ::Val{1}) = Broadcasting(∇)(b.base)
_cp_src(b::CartProdBasis, ::Val{2}) = Broadcasting(∇∇)(b.base)

function _cp_cache(b::CartProdBasis{K}, x, N) where K
  src = _cp_src(b, N)
  cs = return_cache(src, x)
  v = evaluate!(cs, src, x)
  E = representatives_of_componentbasis_dual(VectorValue{K,Float64})
  T = typeof(_cp_insert(first(E), testitem(v)))
  r = CachedArray(zeros(T, size(v, 1), K * size(v, 2)))
  return r, cs, src, E
end

function _cp_eval!(cache, b::CartProdBasis{K}, x, N) where K
  r, cs, src, E = cache
  v = evaluate!(cs, src, x)
  np, nb = size(v)
  setsize!(r, (np, K * nb))
  a = r.array
  @inbounds for c in 1:K
    e = E[c]
    for j in 1:nb, i in 1:np
      a[i, K*(j-1)+c] = _cp_insert(e, v[i, j])
    end
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

# `return_value` must be basis-wise too. Gridap's default route to it is
# `return_value → testargs → testitem → getindex`, i.e. it evaluates one element
# to learn the result type; on a type witness that type is wrong and every cache
# built from it is wrong, silently. Defining it here short-circuits the path.
Gridap.Arrays.return_value(b::CartProdBasis, x::AbstractVector{<:Point}) =
  evaluate(b, x)

Gridap.Arrays.return_value(
  fg::FieldGradientArray{N,<:CartProdBasis}, x::AbstractVector{<:Point}
) where N = evaluate(fg, x)

# Broadcasted differentiation, for the same reason: the generic route asks for
# `testitem`. A basis is atomic, so its derivative is the `FieldGradientArray`.
return_cache(::Broadcasting{typeof(gradient)}, ::CartProdBasis) = nothing
return_cache(::Broadcasting{typeof(∇∇)}, ::CartProdBasis) = nothing

evaluate!(::Nothing, ::Broadcasting{typeof(gradient)}, a::CartProdBasis) =
  FieldGradientArray{1}(a)

evaluate!(::Nothing, ::Broadcasting{typeof(∇∇)}, a::CartProdBasis) =
  FieldGradientArray{2}(a)

############################################################################################
# The stacked DoF basis -- constructed, not wrapped

# `K` copies of the DoF basis `b`, with copy `c` of base DoF `i` at index
# `K*(i-1) + c` and acting on slice `c`:
#
#   σ_{i,c}(u) = σ_i(π_c u),    π_c u = u ⋅ e_c.
#
# The result is an ordinary Gridap DoF basis of the same kind -- the slicing is
# absorbed into the moment weights, which is possible because a moment DoF is
# linear in the field. So nothing here has to differentiate, cache or be
# evaluated field-by-field, and `interpolate` sees a plain `MomentBasedDofBasis`.
function _cp_stack_dofs(b::MomentBasedDofBasis, ::Val{K}) where K
  fm, fom = get_face_moments(b), get_face_own_moments(b)
  E = representatives_of_componentbasis_dual(VectorValue{K,Float64})
  W = typeof(_cp_insert(first(E), zero(eltype(eltype(fm)))))
  moments = map(fm) do M
    ni, nj = size(M)
    A = Array{W}(undef, ni, K * nj)
    for j in 1:nj, c in 1:K, i in 1:ni
      A[i, K*(j-1)+c] = _cp_insert(E[c], M[i, j])
    end
    A
  end
  own = [_cp_expand_dofs(o, K) for o in fom]
  return MomentBasedDofBasis(
    get_nodes(b), moments, get_face_nodes_dofs(b), own, b.operator
  )
end

_cp_stack_dofs(b::ConcatenatedDofVector, ::Val{K}) where K =
  vcat((_cp_stack_dofs(a, Val(K)) for a in b.args)...)

_cp_stack_dofs(b::LinearCombinationDofVector, ::Val{K}) where K =
  linear_combination(_cp_kron(b.values, K), _cp_stack_dofs(b.predofs, Val(K)))

_cp_stack_dofs(b, ::Val) = @notimplemented """\n
Do not know how to stack a $(typeof(b)); add a `_cp_stack_dofs` method for it.
"""

# `kron(C, I_K)` in this file's ordering: entry `C[i,j]` lands at
# `(K*(i-1)+c, K*(j-1)+c)` for every copy `c`.
function _cp_kron!(A::AbstractMatrix, C::AbstractMatrix, K::Integer)
  fill!(A, zero(eltype(A)))
  n, m = size(C)
  @inbounds for j in 1:m, i in 1:n, c in 1:K
    A[K*(i-1)+c, K*(j-1)+c] = C[i, j]
  end
  return A
end

_cp_kron(C::AbstractMatrix, K::Integer) =
  _cp_kron!(zeros(eltype(C), K * size(C, 1), K * size(C, 2)), C, K)

############################################################################################
# The push-forward

# The base push-forward `PF` applied to each of the `K` slices independently,
#
#   φ = Σ_c ι_c( PF(π_c φ̂) ).
#
# This is *not* one of the double Piola maps. Stacking K contravariant Piola maps
# gives φ = det(J)⁻¹ J φ̂, contravariant on the first index only, whereas
# `DoubleContraVariantPiolaMap` is det(J)⁻² J φ̂ Jᵀ. So Hellan-Herrmann-Johnson is
# not obtainable by wrapping Raviart-Thomas, and should not be: the column-wise
# tensor is a different element with a different conformity.
struct CartProdPushforward{K,PF<:Pushforward} <: Pushforward end

_cp_blockwise(pf, v, Jt, E) =
  sum(_cp_insert(e, evaluate(pf, _cp_extract(e, v), Jt)) for e in E)

# `OperationField` builds the cache of its operation once and hands it back to
# every `evaluate!`, so the dual basis is built once per cell, not per value.
return_cache(::CartProdPushforward{K}, ::Number, ::Number) where K =
  representatives_of_componentbasis_dual(VectorValue{K,Float64})

evaluate!(
  E, ::CartProdPushforward{K,PF}, v_ref::Number, Jt::Number
) where {K,PF} = _cp_blockwise(PF(), v_ref, Jt, E)

return_cache(
  ::InversePushforward{CartProdPushforward{K,PF}}, ::Number, ::Number
) where {K,PF} = representatives_of_componentbasis_dual(VectorValue{K,Float64})

evaluate!(
  E, ::InversePushforward{CartProdPushforward{K,PF}}, v_phys::Number, Jt::Number
) where {K,PF} = _cp_blockwise(inverse_map(PF()), v_phys, Jt, E)

############################################################################################
# The reference FE

"""
    struct CartProd{K,N} <: ReferenceFEName end

Reference FE name for `K` stacked copies of an element whose own name type is
`N`. See [`CartProdRefFE`](@ref).
"""
struct CartProd{K,N<:ReferenceFEName} <: ReferenceFEName end

function Pushforward(::Type{CartProd{K,N}}) where {K,N}
  pf = Pushforward(N)
  # K copies of the identity is the identity, and Gridap short-circuits that map
  # rather than evaluating it, so it must be returned unwrapped.
  isa(pf, IdentityPiolaMap) && return IdentityPiolaMap()
  return CartProdPushforward{K,typeof(pf)}()
end

"""
    CartProdRefFE(reffe::ReferenceFE, ::Val{K})
    CartProdRefFE(reffe::ReferenceFE, K::Integer)

`K` independent copies of `reffe`, stacked along a new last index of the value
type — the reference FE of the Cartesian product space `V ⊕ … ⊕ V`, whose values
are

    Float64        ⟹  VectorValue{K}
    VectorValue{d} ⟹  TensorValue{d,K}

so that copy `c` of a vector-valued atom is the `c`-th *column*. A scalar atom
gives vector-valued Morley, Argyris or Hermite; a vector-valued one gives the
column-wise H(div) tensors of elasticity with weakly imposed symmetry [Arnold,
Falk & Winther, Numer. Math. 92 (2002) 401].

The DoFs, shape functions, prebasis and face ownership are all `K` copies of the
base element's, with copy `c` of item `i` at index `K*(i-1) + c`. The conformity
is inherited, and so is the polytope. The base reference FE is kept as the
metadata, which is how `compute_cell_bases_changes` recovers it.

## Examples

    CartProdRefFE(MorleyRefFE(Float64, TRI), Val(2))   # vector-valued Morley
    CartProdRefFE(ReferenceFE(TRI, raviart_thomas, Float64, 1), Val(2))
                                                       # column-wise H(div) tensor
"""
function CartProdRefFE(reffe::ReferenceFE, ::Val{K}) where K
  @notimplementedif K < 1 "Need at least one copy, got $K."
  p = get_polytope(reffe)
  conf = Conformity(reffe)
  nbase = num_dofs(reffe)

  prebasis = CartProdBasis{K}(get_prebasis(reffe))
  shapefuns = CartProdBasis{K}(get_shapefuns(reffe))
  dofs = _cp_stack_dofs(get_dof_basis(reffe), Val(K))

  own = get_face_own_dofs(reffe, conf)
  face_own_dofs = [_cp_expand_dofs(o, K) for o in own]

  N = typeof(get_name(reffe))
  return GenericRefFE{CartProd{K,N}}(
    K * nbase, p, prebasis, dofs, conf, reffe, face_own_dofs, shapefuns
  )
end

CartProdRefFE(reffe::ReferenceFE, K::Integer) = CartProdRefFE(reffe, Val(K))

# The stacked indices of the base DoFs `dofs`: base DoF `d` contributes
# `K*(d-1)+1 … K*d`, kept in base-DoF-major order so that the copy index runs
# fastest, as everywhere here.
_cp_expand_dofs(dofs, K) = Int[K * (d - 1) + c for d in dofs for c in 1:K]

# The base element's permutations, expanded copy-wise: if the base sends its
# local DoF `i` to `p[i]`, the stack sends `K*(i-1)+c` to `K*(p[i]-1)+c`.
# `INVALID_PERM` propagates.
function get_face_own_dofs_permutations(
  reffe::GenericRefFE{CartProd{K,N}}, conf::Conformity
) where {K,N}
  base = get_metadata(reffe)
  base_perms = get_face_own_dofs_permutations(base, conf)
  return [[_cp_expand_perm(p, K) for p in face_perms] for face_perms in base_perms]
end

function _cp_expand_perm(p, K)
  out = fill(INVALID_PERM, K * length(p))
  for i in eachindex(p), c in 1:K
    if p[i] != INVALID_PERM
      out[K*(i-1)+c] = K * (p[i] - 1) + c
    end
  end
  return out
end
