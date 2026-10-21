
struct MaterialHodgeStar{K,D,F<:CellField}
  m::F

  function MaterialHodgeStar{K,D}(m::CellField) where {K,D}
    @check K isa Int && D isa Int """\n
    `MaterialHodgeStar` takes integer `K` and `D`, got $K and $D.
    """
    @check 0 <= K <= D """\n
    A material Hodge star acts on `K`-forms of dimension `D`, which takes
    `0 ≤ K ≤ D`, got K = $K and D = $D.
    """
    x = get_cell_points(get_triangulation(m))
    _check_material_value(Val(K), Val(D), eltype(testitem(m(x))))
    new{K,D,typeof(m)}(m)
  end
end

# The material value must be one the pointwise operator below accepts: the
# `Λ^{D-K}(Cov) ⊗ Λ^K(Con)` tensor itself, or the `L`×`L` matrix of its
# components, whose variance `K` and `D` give.
function _check_material_value(::Val{K}, ::Val{D}, ::Type{Mt}) where {K,D,Mt}
  L = binomial(D,K)
  @check Mt <: MultiValue """\n
  The material of a material Hodge star of `$K`-forms in dimension $D is a
  `MultiValue`, either its $L×$L matrix or the tensor it is read as, got a
  `$Mt`.
  """
  if Mt <: ProductTensorValue
    @check num_indep_components(Mt) == L*L && length(size(Mt)) == D """\n
    The tensor of a material Hodge star of `$K`-forms in dimension $D is a
    `Λ^$(D-K)(Cov) ⊗ Λ^$K(Con)` of $(L*L) independent components, got a `$Mt`.
    """
  else
    sz = size(Mt)
    @check length(sz) == 2 && Int(sz[1]) == L && Int(sz[2]) == L """\n
    The matrix of a material Hodge star of `$K`-forms in dimension $D is
    $L×$L, got a `$Mt` of size $(Int.(sz)).
    """
  end
  nothing
end

_material_tensor_type(::Val{K}, ::Val{D}) where {K,D} =
  ProductTensorValue{ExteriorFormValue{D-K,Cov(D)}, ExteriorVectorValue{K,Con(D)}}

# The matrix and the tensor hold the same components in the same order: the
# independent components of the first factor vary fastest, as a column does.
_as_product(::Val{K}, ::Val{D}, M::ProductTensorValue) where {K,D} = M

_as_product(::Val{K}, ::Val{D}, M::MultiValue{Tuple{L,L},T,2,LL}) where {K,D,L,T,LL} =
  _material_tensor_type(Val(K),Val(D))(M.data)

_as_product(::Val{K}, ::Val{D}, M::MultiValue{S,T,2}) where {K,D,S,T} =
  _material_tensor_type(Val(K),Val(D))(Tuple(get_array(M)))

# ── Constructors ──────────────────────────────────────────────────────────────

# A plain matrix of numbers is one material value, not one value per cell, which
# an `AbstractVector` of values is.
function MaterialHodgeStar{K,D}(M::AbstractMatrix{<:Number}, trian::Triangulation) where {K,D}
  L = binomial(D,K)
  @check size(M) == (L,L) """\n
  The matrix of a material Hodge star of `$K`-forms in dimension $D is $L×$L,
  got $(size(M,1))×$(size(M,2)).
  """
  MaterialHodgeStar{K,D}(TensorValue{L,L}(Tuple(M)), trian)
end

MaterialHodgeStar{K,D}(M, trian::Triangulation) where {K,D} =
  MaterialHodgeStar{K,D}(CellField(M, trian, PhysicalDomain()))

# A product tensor carries `K` and `D`.
function MaterialHodgeStar(
  M::ProductTensorValue{<:ExteriorFormValue{Ko,Dc},<:ExteriorVectorValue{K,Dv}},
  trian::Triangulation) where {Ko,Dc,K,Dv}

  D = Int(Dc)
  @check Int(Dv) == D && Ko == D - K """\n
  A material Hodge star of `$K`-forms in dimension $D returns $(D-K)-forms of
  dimension $D, got a tensor returning $Ko-forms of dimension $(Int(Dc)) from
  $K-vectors of dimension $(Int(Dv)).
  """
  MaterialHodgeStar{K,D}(CellField(M, trian, PhysicalDomain()))
end

"""
    material_hodge_star(M, trian, K, D)
    material_hodge_star(M::CellField, K, D)

The material (constitutive) Hodge star acting on `K`-forms in dimension `D`.
It is a lineal map `Λ^K → Λ^{D-K}`.

A material Hodge star is defined by `L`×`L` components where `L =
binomial(D,K)`, since the map can be represented by a `Λ^{D-K}(Cov) ⊗ Λ^K(Con)`
tensor whose first block indexes the result and whose second block is
contracted with the form argument:

    ⋆ᴹ = material_hodge_star(M, Ω, 1, 2)
    a(u,v) = ∫( v ∧ ⋆ᴹ(u) )dΩ

if `u` is a 2D covector field.

[`hodge_star`](@ref) is the metric Hodge star, such that `M = √|det g|·ε·g⁻¹`
where `ε` is the (anti-symmetric) Levi-Civita symbol.

`M` may be provided as:
- a `TensorValue` or matrix with no variance for a uniform material, or
- an `AbstractVector` of them, for an cell-wise uniform material, or
- a `Function` of the physical coordinate, or a `CellField`.
"""
material_hodge_star(M::CellField, K::Int, D::Int) = MaterialHodgeStar{K,D}(M)
material_hodge_star(M, trian::Triangulation, K::Int, D::Int) = MaterialHodgeStar{K,D}(M, trian)
material_hodge_star(M::ProductTensorValue, trian::Triangulation) = MaterialHodgeStar(M, trian)

# ── Applying the operator ─────────────────────────────────────────────────────

(s::MaterialHodgeStar{K,D})(Mx::MultiValue, ωx::ExteriorFormValue{K}) where {K,D} =
  contracted_product(Val(K), _as_product(Val(K),Val(D),Mx), ωx)

return_cache(::MaterialHodgeStar, ::MultiValue, ::ExteriorFormValue) = nothing
evaluate!(cache, s::MaterialHodgeStar, M::MultiValue, ω::ExteriorFormValue) = s(M, ω)
(s::MaterialHodgeStar)(a::CellField) = Operation(s)(s.m, a)
