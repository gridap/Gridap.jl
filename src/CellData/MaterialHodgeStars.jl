# A material (constitutive) Hodge star is a linear map Λ^K → Λ^{D-K}, the
# composition of the metric Hodge star and an endomorphism of Λ^K
# (pre-composed) or of Λ^{D-K} (post-composed), whose matrix the field `m`
# holds. `inverse` composes with the inverse metric star instead, which is what
# an inverted material Hodge star is made of.
struct MaterialHodgeStar{K,D,PP,F<:CellField}
  m::F
  inverse::Bool

  function MaterialHodgeStar{K,D,PP}(m::CellField, inverse::Bool=false) where {K,D,PP}
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
    new{K,D,PP,typeof(m)}(m, inverse)
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

`⋆ᴹ` and `⋆ᴺ` support `inv`. This inverse is what transform a resistivity to a
conductivity, a magnetic permeability to a magnetic reluctivity, etc.
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
  m = _material_field(M, trian)
  x = get_cell_points(get_triangulation(m))
  _check_material_value(Val(J), Val(D), eltype(testitem(m(x))))
  MaterialHodgeStar{K,D,post_pre}(m)
end

# A material holds the same components at a point whichever domain the field
# carrying them is tagged with, and a physical tag leaves the form the operator
# is applied to in the domain it comes in.
_material_field(M::MultiValue, trian) = CellField(M, trian, PhysicalDomain())
_material_field(M::AbstractArray{<:MultiValue}, trian) =
  CellField(M, trian, PhysicalDomain())
_material_field(M::CellField, trian) = CellField(M, trian, DomainStyle(M))
_material_field(M, trian) = CellField(M, trian)

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

# The degree the material acts in is that of the form `alt_endomorphism` is
# given, `D-K` after the Hodge star and `K` before it.

_metric_star(s::MaterialHodgeStar, a::CellField) =
  s.inverse ? inverse_hodge_star(a) : hodge_star(a)

(s::MaterialHodgeStar{K,D,:post})(a::CellField) where {K,D} =
  Operation(alt_endomorphism)(s.m, _metric_star(s, a))

(s::MaterialHodgeStar{K,D,:pre})(a::CellField) where {K,D} =
  _metric_star(s, Operation(alt_endomorphism)(s.m, a))

# ── Inverting the operator ────────────────────────────────────────────────────

_inverse_material(m::CellField) = Operation(inv)(m)
_inverse_material(m::OperationCellField{DS,typeof(inv)}) where DS = only(m.args)::CellField

# `inv(⋆ᴹ) = M⁻¹ ∘ ⋆⁻¹` and `inv(⋆ᴺ) = ⋆⁻¹ ∘ N⁻¹` map `Λ^{D-K} → Λ^K
inverse_map(s::MaterialHodgeStar{K,D,:pre}) where {K,D} =
  MaterialHodgeStar{D-K,D,:post}(_inverse_material(s.m), !s.inverse)

inverse_map(s::MaterialHodgeStar{K,D,:post}) where {K,D} =
  MaterialHodgeStar{D-K,D,:pre}(_inverse_material(s.m), !s.inverse)

Base.inv(s::MaterialHodgeStar) = inverse_map(s)
