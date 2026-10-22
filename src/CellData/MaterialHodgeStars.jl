# A material (constitutive) Hodge star is a linear map Λ^K → Λ^{D-K}, defined
# by the composition of the metric Hodge star and an endomorphism m/K! on Λ^K
# (pre-composed) or m/(D-K)! on Λ^{D-K} (post-composed).
struct MaterialHodgeStar{K,D,PP,F<:CellField}
  m::F

  function MaterialHodgeStar{K,D,PP}(m::CellField) where {K,D,PP}
    @check K isa Int && D isa Int """\n
    `MaterialHodgeStar` takes integer `K` and `D`, got $K and $D.
    """
    @check 0 <= K <= D """\n
    A material Hodge star acts on `K`-forms of dimension `D`, which takes
    `0 ≤ K ≤ D`, got K = $K and D = $D.
    """
    @check PP === :pre || PP === :post """\n
    The material composes `⋆` either before (`:pre`) or after (`:post`) it, got
    `$PP`.
    """
    new{K,D,PP,typeof(m)}(m)
  end
end

"""
    precomposed_material_hodge_star(M, K, Ω)
    postcomposed_material_hodge_star(N, K, Ω)

A material (or constitutive) Hodge star for the `K`-forms on `Ω`, the linear
map `Λ^K → Λ^{D-K}` resulting from applying `M` before, or `N` after, the
metric Hodge star [`⋆`](@ref):

    ⋆ᴹ = ⋆ ∘ M,    M ∈ End(Λ^K),
    ⋆ᴺ = N ∘ ⋆,    N ∈ End(Λ^{D-K}).

# Arguments

`D` is the manifold dimension of a triangulation `Ω`.

`M`/`N` is the `L`×`L` matrix or tensor of the corresponding endomorphism, with
the components ordered similarly to those of `K` or `D-K` form, respectively.
These act as

    (M(ω))ᵢ = Σⱼ Mᵢⱼ ωⱼ,    (N(ω))ᵢ = Σⱼ Nᵢⱼ ωⱼ.

`M`/`N` may be given as:
- a `MultiValue` or an `AbstractMatrix` with no variance, for a uniform material,
- an `AbstractVector` of them, for a cell-wise uniform material,
- a `Function` of the physical coordinate, or a `CellField`.

# Extended help

### Relations with the inner product of forms

The definition of the [`hodge_star`](@ref) in terms of the inner product of
`K`-forms `ω` and `η`

    ω ∧ ⋆η = (ω ⨟ η) μ,  μ = √|det g| dx¹∧…∧dx^D

generalises to

    ω ∧ ⋆ᴹ(η) =  (  ω ⨟ M(η)  ) μ,
    ω ∧ ⋆ᴺ(η) = s( ⋆ω ⨟ N(⋆η) ) μ,

where `s = sgn(det(g))`. ⋆ᴹ` and `⋆ᴺ` are the same operator if `M` and `N` verify

    M = ⋆⁻¹ ∘ N ∘ ⋆,     N = ⋆ ∘ M ∘ ⋆⁻¹.

### Inverse

A `⋆` squares to `s(-1)^(K*(D-K))` on `Λ^K`, so inverting a material
Hodge star writes

    inv(⋆ᴹ) = M⁻¹ ∘ ⋆⁻¹ = s(-1)^(K*(D-K)) M⁻¹ ∘ ⋆,
    inv(⋆ᴺ) = ⋆⁻¹ ∘ N⁻¹ = s(-1)^(K*(D-K)) ⋆ ∘ N⁻¹.

This inverse is what transform a resistivity to a conductivity, a magnetic
permeability to a magnetic reluctivity, etc.
"""
postcomposed_material_hodge_star(N, K::Int, Ω::Triangulation) =
  _material_hodge_star(N, K, Ω, :post)

"""
    precomposed_material_hodge_star(M, K, Ω)

See [`postcomposed_material_hodge_star`](@ref)
"""
precomposed_material_hodge_star(M, K::Int, Ω::Triangulation) =
  _material_hodge_star(M, K, Ω, :pre)

_material_degree(K::Int, D::Int, post_pre::Symbol) = post_pre === :post ? D-K : K

# A plain matrix of numbers is one material value, not one value per cell, which
# an `AbstractVector` of values is.
function _material_hodge_star(M::AbstractMatrix{<:Number}, K, trian, post_pre)
  D = num_cell_dims(trian)
  L = binomial(D, _material_degree(K,D,post_pre))
  @check size(M) == (L,L) """\n
  The material of this material Hodge star is $L×$L, got $(size(M,1))×$(size(M,2)).
  """
  _material_hodge_star(TensorValue{L,L}(Tuple(M)), K, trian, post_pre)
end

function _material_hodge_star(M, K, trian, post_pre)
  D = num_cell_dims(trian)
  @check 0 <= K <= D """\n
  A material Hodge star acts on `K`-forms of dimension `D`, which takes
  `0 ≤ K ≤ D`, got K = $K and D = $(D).
  """
  J = _material_degree(K,D,post_pre)
  m = CellField(M, trian)
  x = get_cell_points(get_triangulation(m))
  _check_material_value(Val(J), Val(D), eltype(testitem(m(x))))
  MaterialHodgeStar{K,D,post_pre}(_scaled_material(M, trian, Val(factorial(J))))
end

# The `J!` that the stored component contraction of a `J`-form carries is
# divided out of the material once, here, so that applying the operator
# contracts the material components as they are stored.
_scaled(::Val{f}, M) where f = M/f

_scaled_material(M::MultiValue, trian, f::Val) = CellField(_scaled(f,M), trian)
_scaled_material(M::AbstractArray{<:MultiValue}, trian, f::Val) =
  CellField(map(m -> _scaled(f,m), M), trian)
_scaled_material(M::Function, trian, f::Val) = CellField(x -> _scaled(f,M(x)), trian)
_scaled_material(M, trian, f::Val) = Operation(m -> _scaled(f,m))(CellField(M, trian))

# The material is the `L`×`L` matrix of an endomorphism of `Λ^J`, nothing else:
# the `Λ^{D-K}(Cov) ⊗ Λ^K(Con)` tensor of the whole map is built from it.
function _check_material_value(::Val{J}, ::Val{D}, ::Type{Mt}) where {J,D,Mt}
  L = binomial(D,J)
  @check Mt <: MultiValue """\n
  The material of a material Hodge star is the $L×$L matrix of an endomorphism
  of the $J-forms of dimension $D, a `MultiValue`, got a `$Mt`.
  """
  sz = size(Mt)
  @check length(sz) == 2 && Int(sz[1]) == L && Int(sz[2]) == L """\n
  The material of a material Hodge star is the $L×$L matrix of an endomorphism
  of the $J-forms of dimension $D, got a `$Mt` of size $(Int.(sz)).
  """
  nothing
end

# ── Applying the operator ─────────────────────────────────────────────────────

# Application of the endomorphism `m` to the J-form value, is a J contraction of `m/J!`.
# The contraction is optimized by a specific method on ProductTensorValue

function _as_endomorphism(
  ::Val{J}, ::Val{D}, M::MultiValue{Tuple{L,L}}) where {J,D,L}

  ET = ProductTensorValue{ExteriorFormValue{J,Cov(D)}, ExteriorVectorValue{J,Con(D)}}
  ET(Tuple(get_array(M)))
end

(::MaterialHodgeStar{K,D,:post})(Mx::MultiValue, ωx::ExteriorFormValue) where {K,D} =
  contracted_product(Val(D-K), _as_endomorphism(Val(D-K),Val(D),Mx), ωx)

(::MaterialHodgeStar{K,D,:pre})(Mx::MultiValue, ωx::ExteriorFormValue) where {K,D} =
  contracted_product(Val(K), _as_endomorphism(Val(K),Val(D),Mx), ωx)

return_cache(::MaterialHodgeStar, ::MultiValue, ::ExteriorFormValue) = nothing
evaluate!(cache, s::MaterialHodgeStar, M::MultiValue, ω::ExteriorFormValue) = s(M, ω)

(s::MaterialHodgeStar{K,D,:post})(a::CellField) where {K,D} =
  Operation(s)(s.m, hodge_star(a))

(s::MaterialHodgeStar{K,D,:pre})(a::CellField) where {K,D} =
  hodge_star(Operation(s)(s.m, a))
