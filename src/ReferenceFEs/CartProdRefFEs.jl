# The Cartesian product V₁ × … × V_K of reference elements; see `CartProdRefFE`.
#
# We stack the factors along the LAST index. Two reasons:
#   - Derivatives in Gridap are ∇A[k,i,j] = ∂_k A_ij. By stacking on the last index, 
#     the derivative and the stacking (injection operator below) commute.
#   - The last index is also the slow index in Julia. This means that using this 
#     convention, each copy of the factor is stored contiguously in memory. Thus, 
#     block access patterns are more efficient.
# 
# The injection ι_c(v) = v ⊗ e_c and its inverse π_c(v) = v ⋅ e_c are the two
# constant linear maps the whole construction rests on:
@inline _cp_insert(e, v) = outer(v, e)
@inline _cp_extract(e, v) = v ⋅ e
# where `v` is the element of `V` to be injected or extracted 
# and `e` is an element of the orthonormal basis of the stacked product space.
# If `e` is the `c`-th element of the orthonormal basis, these operators 
# inject into and extract from the `c`-th block.

############################################################################################
# Shape functions 
# 
# For each component `c` and each basis function `i` in the `c`-th factor basis, we have
#  
# Φ_{i,c}(x) = ι_c(Φ_i(x))
# 
# and since derivatives commute with the injection, we also have
#
# ∂_k Φ_{i,c}(x) = ι_c(∂_k Φ_i(x))
#

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

_cp_src(f, ::Val{0}) = f
_cp_src(f, ::Val{1}) = Broadcasting(∇)(f)
_cp_src(f, ::Val{2}) = Broadcasting(∇∇)(f)

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
# DoF basis
#
# For each component `c` and each basis function `i` in the `c`-th factor dof basis, we have
#
#   σ_{i,c}(u) = σ_i(π_c u),    π_c u = u ⋅ e_c,
#

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

# ... except for a nodal basis, which stacks whole with no slicing at all. 
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
#
# Pushforwards apply component-wise, that is for each component `c` and each 
# index `i` of the `c`-th component, we have
#
#   Φ_c = ι_c(PF_c(π_c(ϕ))), that is Φ_{c,i} = ι_c(PF_c(ϕ_{c,i}))
#

struct CartProdPushforward{K,PFS<:Tuple} <: Pushforward end

# `IdentityPiolaMap` is short-circuited for efficiency.
@inline _cp_apply(::IdentityPiolaMap, v, Jt) = v
@inline _cp_apply(::InversePushforward{IdentityPiolaMap}, v, Jt) = v
@inline _cp_apply(pf, v, Jt) = evaluate(pf, v, Jt)

_cp_maps(::Type{PFS}) where PFS<:Tuple =
  ntuple(i -> fieldtype(PFS, i)(), fieldcount(PFS))

_cp_blockwise(pfs, v, Jt, E) =
  sum(_cp_insert(E[i], _cp_apply(pfs[i], _cp_extract(E[i], v), Jt)) for i in eachindex(E))

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

Cartesian product of reference elements. We provide two constructors:

-  `CartProdRefFE(reffe::ReferenceFE, K::Integer)` builds `V^K`, the Cartesian 
   product of `K` copies of a single space `V`, given by `reffe`.

-  `CartProdRefFE(reffes::ReferenceFE...)` builds `V_1 × … × V_K`, the Cartesian
    product of `K` different spaces, given by the `reffes` tuple.

By design, the individual spaces are stacked along the last index of the new value type.
This is deliberately chosen to be compatible with Gridap's conventions for derivatives,
and to keep DoFs of each factor in contiguous blocks.

The factors must share a polytope and a value type.

## Examples:

    CartProdRefFE(MorleyRefFE(Float64, TRI), 2)        # vector-valued Morley
    CartProdRefFE(ReferenceFE(TRI, raviart_thomas, Float64, 1), 2)
                                                       # column-wise H(div) tensor
    CartProdRefFE(ReferenceFE(TRI, lagrangian, Float64, 2),
                  ReferenceFE(TRI, lagrangian, Float64, 2),
                  ArgyrisRefFE(Float64, TRI))          # Kirchhoff--Love shell

"""
struct CartProd{K,NS} <: ReferenceFEName end

"""
    struct CartProdConformity{K} <: Conformity

    CartProdConformity(confs::Conformity...)
    CartProdConformity(conf::Conformity, K::Integer)

The conformity of a [`CartProd`](@ref) element whose factors do not all share one.
"""
struct CartProdConformity{K,C<:NTuple{K,Conformity}} <: Conformity
  confs::C
end

CartProdConformity(confs::Conformity...) = CartProdConformity(confs)
CartProdConformity(conf::Conformity, K::Integer) =
  CartProdConformity(ntuple(i -> conf, K))

# only what every factor accepts, which for a mixed product is just `:L2`
valid_conformity_symbols(conf::CartProdConformity) =
  Tuple(intersect(map(valid_conformity_symbols, conf.confs)...))

_cp_factor_names(::Val{K}, ::Type{N}) where {K,N<:ReferenceFEName} = ntuple(i -> N, K)
_cp_factor_names(::Val{K}, ::Type{NS}) where {K,NS<:Tuple} = ntuple(i -> fieldtype(NS, i), K)

function Pushforward(N::Type{CartProd{K,NS}}, conf::CartProdConformity{K}) where {K,NS}
  pfs = map(Pushforward, _cp_factor_names(Val(K), NS), conf.confs)
  # K copies of the identity is the identity, and Gridap short-circuits that map
  # rather than evaluating it, so it must be returned unwrapped.
  all(pf -> isa(pf, IdentityPiolaMap), pfs) && return IdentityPiolaMap()
  return CartProdPushforward{K,Tuple{map(typeof, pfs)...}}()
end

Pushforward(N::Type{CartProd{K,NS}}, conf::Conformity) where {K,NS} =
  Pushforward(N, CartProdConformity(conf, K))

# needed only to break the tie with the generic `L2Conformity` method
Pushforward(N::Type{CartProd{K,NS}}, conf::L2Conformity) where {K,NS} =
  Pushforward(N, CartProdConformity(conf, K))

"""

    CartProdRefFE(reffes::ReferenceFE...)

Given `K` local spaces given by the reference finite elements `reffes`, 
returns the Cartesian product reference finite element `V_1 × … × V_K`.
See [`CartProd`](@ref) for details.

"""
function CartProdRefFE(reffes::ReferenceFE...)
  K = length(reffes)
  @notimplementedif K < 1 "Need at least one factor."
  @notimplementedif any(r -> get_polytope(r) != get_polytope(first(reffes)), reffes) """\n
  All factors must share a polytope.
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

"""
    CartProdRefFE(reffe::ReferenceFE, K::Integer)

Given `V` a local space given by the reference finite element `reffe`, and `K` a positive integer,
returns the Cartesian product reference finite element `V^K = V × … × V`.
See [`CartProd`](@ref) for details.

"""
CartProdRefFE(reffe::ReferenceFE, K::Integer) = CartProdRefFE(reffe, Val(K))

_cp_value_type(reffe) = typeof(testitem(
  evaluate(get_shapefuns(reffe), get_vertex_coordinates(get_polytope(reffe)))
))

# The factors never mix, so everything is blocked: the DoFs are the `vcat` of the
# slots, and a face owns each factor's own DoFs shifted past the factors before
# it.
#
# The conformity is the factors' when they all share one, so that a homogeneous
# product looks like its factors to everything that dispatches on it, and a
# `CartProdConformity` otherwise.
function _cp_reffe(::Type{Name}, reffes, prebasis, shapefuns) where Name
  K = length(reffes)
  p = get_polytope(first(reffes))
  confs = map(Conformity, reffes)
  conf = all(==(first(confs)), confs) ? first(confs) : CartProdConformity(confs)

  dofs = _cp_stack_dofs(map(get_dof_basis, reffes), Val(K))
  face_own_dofs = _cp_face_own_dofs(reffes, CartProdConformity(confs))

  return GenericRefFE{Name}(
    sum(num_dofs, reffes), p, prebasis, dofs, conf, reffes, face_own_dofs, shapefuns
  )
end

# `_cp_offsets(ns)[c]` is the number of items before factor `c`
_cp_offsets(ns) = (0, cumsum(ns)[1:end-1]...)

# Ownership and permutations are the factors', each under its own conformity,
# blocked. Any other conformity is applied to every factor alike, which is how
# `conformity=:L2` detaches a product.

function _cp_face_own_dofs(reffes, conf::CartProdConformity)
  offs = _cp_offsets(map(num_dofs, reffes))
  owns = map(get_face_own_dofs, reffes, conf.confs)
  [Int[offs[c] + d for c in eachindex(owns) for d in owns[c][f]]
   for f in eachindex(first(owns))]
end

get_face_own_dofs(
  reffe::GenericRefFE{CartProd{K,NS}}, conf::CartProdConformity{K}
) where {K,NS} = _cp_face_own_dofs(get_metadata(reffe), conf)

get_face_own_dofs(reffe::GenericRefFE{CartProd{K,NS}}, conf::Conformity) where {K,NS} =
  get_face_own_dofs(reffe, CartProdConformity(conf, K))

# The offsets are per face: factor `c`'s own DoFs on a face sit after those of
# the factors before it *on that same face*. `INVALID_PERM` propagates.
function get_face_own_dofs_permutations(
  reffe::GenericRefFE{CartProd{K,NS}}, conf::CartProdConformity{K}
) where {K,NS}
  reffes = get_metadata(reffe)
  perms = map(get_face_own_dofs_permutations, reffes, conf.confs)
  owns = map(get_face_own_dofs, reffes, conf.confs)
  map(eachindex(first(perms))) do f
    offs = _cp_offsets(map(o -> length(o[f]), owns))
    map(eachindex(first(perms)[f])) do pindex
      Int[q == INVALID_PERM ? INVALID_PERM : offs[c] + q
          for c in 1:K for q in perms[c][f][pindex]]
    end
  end
end

get_face_own_dofs_permutations(
  reffe::GenericRefFE{CartProd{K,NS}}, conf::Conformity
) where {K,NS} = get_face_own_dofs_permutations(reffe, CartProdConformity(conf, K))

################################################################################
# Change of basis
#
# The cell-local map. The mesh-level `compute_cell_bases_changes`, which reads
# the orientation data off the model and maps this over the cells, lives in
# src/FESpaces/Pullbacks.jl.

# The factors of a stacked element never mix, so its change of basis is the
# factors' assembled block diagonally, in the blocked ordering of
# `CartProdRefFEs.jl`. A factor that needs no change of basis contributes an
# identity block; when none of them do, `nothing` propagates.

struct CartProdBlockDiag <: Map end

function return_cache(::CartProdBlockDiag, Ps::AbstractMatrix...)
  T = promote_type(map(eltype, Ps)...)
  CachedArray(zeros(T, sum(P -> size(P, 1), Ps), sum(P -> size(P, 2), Ps)))
end

function evaluate!(cache, ::CartProdBlockDiag, Ps::AbstractMatrix...)
  setsize!(cache, (sum(P -> size(P, 1), Ps), sum(P -> size(P, 2), Ps)))
  A = cache.array
  fill!(A, zero(eltype(A)))
  io = jo = 0
  for P in Ps
    ni, mj = size(P)
    @inbounds A[io+1:io+ni, jo+1:jo+mj] .= P
    io += ni
    jo += mj
  end
  return A
end

_cp_maps(::CartProdPushforward{K,PFS}, ::Val{K}) where {K,PFS} =
  _cp_maps(PFS)

_cp_maps(pf::IdentityPiolaMap, ::Val{K}) where K = ntuple(i -> pf, K)

################################################################################
# DOF scaling
#
# The factors never mix, so each factor's own setter scales its block through a
# view.
function get_dofscale_setter_function(
  reffe::GenericRefFE{CartProd{K,NS}}, pf::Pushforward
) where {K,NS}
  reffes = get_metadata(reffe)
  ndofs = map(num_dofs, reffes)
  blocks = map((off, n) -> off .+ (1:n), _cp_offsets(ndofs), ndofs)
  setters = map(get_dofscale_setter_function, reffes, _cp_maps(pf, Val(K)))
  factor_own_dofs = map(get_face_own_dofs, reffes)
  blocked_own_dofs = get_face_own_dofs(reffe)

  let setters=setters, blocks=blocks, factor_own_dofs=factor_own_dofs,
      blocked_own_dofs=blocked_own_dofs
    @inline function(dofscale, face_own_dofs, face_meshsize)
      @check face_own_dofs == blocked_own_dofs "unexpected face ownership in the DoF scale setter. scale_dof may be disabled as a temporary solution."
      foreach(setters, blocks, factor_own_dofs) do setter, block, own_dofs
        setter(view(dofscale, block), own_dofs, face_meshsize)
      end
      nothing
    end
  end
end

