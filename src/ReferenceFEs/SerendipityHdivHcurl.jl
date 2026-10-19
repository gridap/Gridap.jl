################################################################################
# Serendipity-type H(div) and H(curl) elements on quadrilaterals and hexahedra:
# Arnold-Awanou's 𝓢ᵣΛᵏ and Gillette-Kloefkorn's trimmed 𝓢⁻ᵣΛᵏ, for k = 1, D-1.
#
# References:
# [AA] D.N. Arnold & G. Awanou, Finite element differential forms on cubical
#      meshes, Math. Comp. 83 (2014) 1551-1570.
# [GK] A. Gillette & T. Kloefkorn, Trimmed serendipity finite element
#      differential forms, Math. Comp. 88 (2019) 583-606.

"""
    struct SerendipityDiv <: ReferenceFEName end

Reference FE name for the Arnold--Awanou H(div) element on n-cubes. See
[`SerendipityDivRefFE`](@ref). Its singleton instance is [`sdiv`](@ref).
"""
struct SerendipityDiv <: ReferenceFEName end

"""
    const sdiv = SerendipityDiv()

Singleton of the [`SerendipityDiv`](@ref) reference FE name.
"""
const sdiv = SerendipityDiv()

"""
    struct SerendipityCurl <: ReferenceFEName end

Reference FE name for the Arnold--Awanou H(curl) element on n-cubes. See
[`SerendipityCurlRefFE`](@ref). Its singleton instance is [`scurl`](@ref).
"""
struct SerendipityCurl <: ReferenceFEName end

"""
    const scurl = SerendipityCurl()

Singleton of the [`SerendipityCurl`](@ref) reference FE name.
"""
const scurl = SerendipityCurl()

"""
    struct TrimmedSerendipityDiv <: ReferenceFEName end

Reference FE name for the Gillette--Kloefkorn trimmed serendipity H(div)
element on n-cubes. See [`TrimmedSerendipityDivRefFE`](@ref). Its singleton
instance is [`trimmed_sdiv`](@ref).
"""
struct TrimmedSerendipityDiv <: ReferenceFEName end

"""
    const trimmed_sdiv = TrimmedSerendipityDiv()

Singleton of the [`TrimmedSerendipityDiv`](@ref) reference FE name.
"""
const trimmed_sdiv = TrimmedSerendipityDiv()

"""
    struct TrimmedSerendipityCurl <: ReferenceFEName end

Reference FE name for the Gillette--Kloefkorn trimmed serendipity H(curl)
element on n-cubes. See [`TrimmedSerendipityCurlRefFE`](@ref). Its singleton
instance is [`trimmed_scurl`](@ref).
"""
struct TrimmedSerendipityCurl <: ReferenceFEName end

"""
    const trimmed_scurl = TrimmedSerendipityCurl()

Singleton of the [`TrimmedSerendipityCurl`](@ref) reference FE name.
"""
const trimmed_scurl = TrimmedSerendipityCurl()

const _SerendipityHdivHcurl = Union{
  SerendipityDiv,SerendipityCurl,TrimmedSerendipityDiv,TrimmedSerendipityCurl}

Pushforward(::Type{SerendipityDiv}) = ContraVariantPiolaMap()
Pushforward(::Type{SerendipityCurl}) = CoVariantPiolaMap()
Pushforward(::Type{TrimmedSerendipityDiv}) = ContraVariantPiolaMap()
Pushforward(::Type{TrimmedSerendipityCurl}) = CoVariantPiolaMap()

"""
    SerendipityDivRefFE(::Type{T}, p::Polytope{D}, r::Integer)

The serendipity H(div) reference FE of degree `r ≥ 1` on the `D`-cube `p`, with
`T` the scalar type, for `D = 2` and `D = 3`:

    𝓢ᵣΛᴰ⁻¹(□ᴰ) = 𝓟ᵣ(□ᴰ;Rᴰ) ⊕ d𝓙ᵣ₊₁Λᴰ⁻²

of dimension `r² + 3r + 4` on a square and `(r+1)(r² + 5r + 12)/2` on a cube
[AA, (17) and Thm. 3.6]. On a square it is the rectangular BDM element [Brezzi,
Douglas & Marini, Numer. Math. 47 (1985) 217]. The element is div-conforming;
`div` maps it onto `𝓟ᵣ₋₁(□ᴰ)`.

# Extended help

## Prebasis

`FEEC_poly_basis(Val(D),T,r,D-1,:S)`, see [`SerendipityPolyBasis`](@ref).

## Moments

With `n` the unit outward normal of the facet `f`,

    ℓ^{f,q}(u) = ∫_f (u⋅n) q,     q ∈ 𝓟ᵣ(f),            for each facet f
    ℓ^{p}(u)   = ∫_K u⋅p,         p ∈ 𝓟ᵣ₋₂(□ᴰ;Rᴰ),

that is `2D dim 𝓟ᵣ(f)` facet DoFs plus `D dim 𝓟ᵣ₋₂(□ᴰ)` interior ones.
"""
function SerendipityDivRefFE(::Type{T}, p::Polytope{D}, r::Integer) where {T,D}
  @notimplementedif !(D in (2,3) && is_n_cube(p)) """\n
  The serendipity H(div) element is only implemented on squares and cubes, got $p.
  """
  @notimplementedif r < 1 "The serendipity H(div) element needs r ≥ 1, got $r."

  prebasis = FEEC_poly_basis(Val(D), T, r, D-1, :S; rotate_90=(D==2))

  function cmom(φ,μ,ds) # Cell moment function: σ_K(φ,μ) = ∫(φ·μ)dK
    Broadcasting(Operation(⋅))(φ,μ)
  end
  function fmom(φ,μ,ds) # Face moment function: σ_F(φ,μ) = ∫((φ·n)*μ)dF, n outward
    n = get_facet_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ,n)
    Broadcasting(Operation(*))(φn,μ)
  end

  fb = LegendreBasis(Val(D-1), T, r, Polynomials._p_filter)  # 𝓟ᵣ(f)
  moments = Tuple[(get_dimrange(p, D-1), fmom, fb)]
  if r > 1
    cb = LegendreBasis(Val(D), VectorValue{D,T}, r-2, Polynomials._p_filter) # 𝓟ᵣ₋₂(K;Rᴰ)
    push!(moments, (get_dimrange(p, D), cmom, cb))
  end
  conf = DivConformity()
  dofs = MomentBasedDofBasis(p, prebasis, moments)
  face_own_dofs_permutations = moment_face_own_dofs_permutations(p, moments, conf)
  GenericRefFE{SerendipityDiv}(
    length(dofs), p, prebasis, dofs, conf, nothing, get_face_own_moments(dofs);
    face_own_dofs_permutations)
end

function ReferenceFE(p::Polytope, ::SerendipityDiv, ::Type{T}, order) where T
  SerendipityDivRefFE(T, p, order)
end

"""
    SerendipityCurlRefFE(::Type{T}, p::Polytope{D}, r::Integer)

The serendipity H(curl) reference FE of degree `r ≥ 1` on the `D`-cube `p`, with
`T` the scalar type, for `D = 2` and `D = 3`:

    𝓢ᵣΛ¹(□ᴰ) = 𝓟ᵣ(□ᴰ;Rᴰ) ⊕ 𝓙ᵣΛ¹ ⊕ d𝓙ᵣ₊₁Λ⁰

of dimension `Σ_{d=1}^{min(D,⌊r/2⌋+1)} 2ᴰ⁻ᵈ C(D,d) C(r-d+2,d) d` [AA, (17) and
(22)], `r² + 3r + 4` on a square. The element is curl-conforming.

# Extended help

## Prebasis

`FEEC_poly_basis(Val(D),T,r,1,:S)`, see [`SerendipityPolyBasis`](@ref).

## Moments

With `t` the unit tangent of the edge `e` and `n` the unit outward normal of the
face `f`,

    ℓ^{e,q}(u) = ∫_e (u⋅t) q,     q ∈ 𝓟ᵣ(e),                for each edge e
    ℓ^{f,p}(u) = ∫_f (u×n)⋅p,     p ∈ 𝓟ᵣ₋₂(f;R²),           for each face f (3D)
    ℓ^{p}(u)   = ∫_K u⋅p,         p ∈ 𝓟ᵣ₋₂₍ᴰ₋₁₎(□ᴰ;Rᴰ).
"""
function SerendipityCurlRefFE(::Type{T}, p::Polytope{D}, r::Integer) where {T,D}
  @notimplementedif !(D in (2,3) && is_n_cube(p)) """\n
  The serendipity H(curl) element is only implemented on squares and cubes, got $p.
  """
  @notimplementedif r < 1 "The serendipity H(curl) element needs r ≥ 1, got $r."

  prebasis = FEEC_poly_basis(Val(D), T, r, 1, :S)

  function cmom(φ,μ,ds) # Cell moment function: σ_K(φ,μ) = ∫(φ·μ)dK
    Broadcasting(Operation(⋅))(φ,μ)
  end
  function fmom(φ,μ,ds) # Face moment function: σ_F(φ,μ) = ∫((φ×n)·Eμ)dF, n outward
    n = get_facet_normal(ds)
    E = get_extension(ds)
    Eμ = Broadcasting(Operation(⋅))(E,μ)
    φn = Broadcasting(Operation(×))(φ,n)
    Broadcasting(Operation(⋅))(φn,Eμ)
  end
  function emom(φ,μ,ds) # Edge moment function: σ_E(φ,μ) = ∫((φ·t)*μ)dE
    t = get_edge_tangent(ds)
    φt = Broadcasting(Operation(⋅))(φ,t)
    Broadcasting(Operation(*))(φt,μ)
  end

  eb = LegendreBasis(Val(1), T, r, Polynomials._p_filter)  # 𝓟ᵣ(e)
  moments = Tuple[(get_dimrange(p, 1), emom, eb)]
  if D == 3 && r > 1
    fb = LegendreBasis(Val(2), VectorValue{2,T}, r-2, Polynomials._p_filter) # 𝓟ᵣ₋₂(f;R²)
    push!(moments, (get_dimrange(p, 2), fmom, fb))
  end
  ρ = r - 2*(D-1)
  if ρ ≥ 0
    cb = LegendreBasis(Val(D), VectorValue{D,T}, ρ, Polynomials._p_filter)   # 𝓟ρ(K;Rᴰ)
    push!(moments, (get_dimrange(p, D), cmom, cb))
  end
  conf = CurlConformity()
  dofs = MomentBasedDofBasis(p, prebasis, moments)
  face_own_dofs_permutations = moment_face_own_dofs_permutations(p, moments, conf)
  GenericRefFE{SerendipityCurl}(
    length(dofs), p, prebasis, dofs, conf, nothing, get_face_own_moments(dofs);
    face_own_dofs_permutations)
end

function ReferenceFE(p::Polytope, ::SerendipityCurl, ::Type{T}, order) where T
  SerendipityCurlRefFE(T, p, order)
end

"""
    TrimmedSerendipityDivRefFE(::Type{T}, p::Polytope{D}, k::Integer)

The trimmed serendipity H(div) reference FE of order `k ≥ 0` on the `D`-cube
`p`, with `T` the scalar type, for `D = 2` and `D = 3`:

    𝓢⁻ₖ₊₁Λᴰ⁻¹(□ᴰ) = 𝓢ₖΛᴰ⁻¹ + κ𝓢ₖΛᴰ

[GK, (3.1) and Table 1]. It contains `𝓟ₖ(□ᴰ;Rᴰ)` and is the lowest order
Raviart-Thomas element at `k = 0`. The element is div-conforming; `div` maps
it onto `𝓟ₖ(□ᴰ)`.

# Extended help

## Prebasis

`FEEC_poly_basis(Val(D),T,k+1,D-1,:S⁻)`, see [`SerendipityPolyBasis`](@ref).

## Moments

With `n` the unit outward normal of the facet `f`,

    ℓ^{f,q}(u) = ∫_f (u⋅n) q,     q ∈ 𝓟ₖ(f),                    for each facet f
    ℓ^{p}(u)   = ∫_K u⋅p,         p ∈ 𝓟ₖ₋₂(□ᴰ;Rᴰ) ⊕ ∇𝓗ₖ(□ᴰ),

[GK, (4.1)].
"""
function TrimmedSerendipityDivRefFE(::Type{T}, p::Polytope{D}, k::Integer) where {T,D}
  @notimplementedif !(D in (2,3) && is_n_cube(p)) """\n
  The trimmed serendipity H(div) element is only implemented on squares and cubes, got $p.
  """
  @notimplementedif k < 0 "The trimmed serendipity H(div) element needs k ≥ 0, got $k."

  prebasis = FEEC_poly_basis(Val(D), T, k+1, D-1, :S⁻; rotate_90=(D==2))

  function cmom(φ,μ,ds) # Cell moment function: σ_K(φ,μ) = ∫(φ·μ)dK
    Broadcasting(Operation(⋅))(φ,μ)
  end
  function fmom(φ,μ,ds) # Face moment function: σ_F(φ,μ) = ∫((φ·n)*μ)dF, n outward
    n = get_facet_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ,n)
    Broadcasting(Operation(*))(φn,μ)
  end

  fb = LegendreBasis(Val(D-1), T, k, Polynomials._p_filter)  # 𝓟ₖ(f)
  moments = Tuple[(get_dimrange(p, D-1), fmom, fb)]
  if k > 1
    cb = LegendreBasis(Val(D), VectorValue{D,T}, k-2, Polynomials._p_filter) # 𝓟ₖ₋₂(K;Rᴰ)
    push!(moments, (get_dimrange(p, D), cmom, cb))
  end
  if k > 0
    gb = Polynomials._exact_poly_basis(Legendre, T, D, k-1, 1)               # ∇𝓗ₖ(K)
    push!(moments, (get_dimrange(p, D), cmom, gb))
  end
  conf = DivConformity()
  dofs = MomentBasedDofBasis(p, prebasis, moments)
  face_own_dofs_permutations = moment_face_own_dofs_permutations(p, moments, conf)
  GenericRefFE{TrimmedSerendipityDiv}(
    length(dofs), p, prebasis, dofs, conf, nothing, get_face_own_moments(dofs);
    face_own_dofs_permutations)
end

function ReferenceFE(p::Polytope, ::TrimmedSerendipityDiv, ::Type{T}, order) where T
  TrimmedSerendipityDivRefFE(T, p, order)
end

"""
    TrimmedSerendipityCurlRefFE(::Type{T}, p::Polytope{D}, k::Integer)

The trimmed serendipity H(curl) reference FE of order `k ≥ 0` on the `D`-cube
`p`, with `T` the scalar type, for `D = 2` and `D = 3`:

    𝓢⁻ₖ₊₁Λ¹(□ᴰ) = 𝓢ₖΛ¹ + κ𝓢ₖΛ²

[GK, (3.1) and Table 1]. It contains `𝓟ₖ(□ᴰ;Rᴰ)` and is the lowest order
Nédélec element at `k = 0`. The element is curl-conforming.

# Extended help

## Prebasis

`FEEC_poly_basis(Val(D),T,k+1,1,:S⁻)`, see [`SerendipityPolyBasis`](@ref).

## Moments

With `t` the unit tangent of the edge `e` and `n` the unit outward normal of the
face `f`,

    ℓ^{e,q}(u) = ∫_e (u⋅t) q,   q ∈ 𝓟ₖ(e),                       for each edge e
    ℓ^{f,p}(u) = ∫_f (u×n)⋅p,   p ∈ 𝓟ₖ₋₂(f;R²) ⊕ ∇𝓗ₖ(f),          for each face f (3D)
    ℓ^{p}(u)   = ∫_K u⋅p,       p ∈ 𝓟ₖ₋₂(□²;R²) ⊕ curl 𝓗ₖ(□²)          (2D)
                                p ∈ 𝓟ₖ₋₄(□³;R³) ⊕ curl 𝓗ₖ₋₂(□³;R³)     (3D)

[GK, (4.1)]: on a `d`-face, the index space `𝓟ₛ₋₁Λᵈ⁻¹ ⊕ d𝓗ₛ₊₁Λᵈ⁻²`,
`s = k+1-2(d-1)`, paired with the trace of `u` by the wedge product.
"""
function TrimmedSerendipityCurlRefFE(::Type{T}, p::Polytope{D}, k::Integer) where {T,D}
  @notimplementedif !(D in (2,3) && is_n_cube(p)) """\n
  The trimmed serendipity H(curl) element is only implemented on squares and cubes, got $p.
  """
  @notimplementedif k < 0 "The trimmed serendipity H(curl) element needs k ≥ 0, got $k."

  prebasis = FEEC_poly_basis(Val(D), T, k+1, 1, :S⁻)

  function cmom(φ,μ,ds) # Cell moment function: σ_K(φ,μ) = ∫(φ·μ)dK
    Broadcasting(Operation(⋅))(φ,μ)
  end
  function fmom(φ,μ,ds) # Face moment function: σ_F(φ,μ) = ∫((φ×n)·Eμ)dF, n outward
    n = get_facet_normal(ds)
    E = get_extension(ds)
    Eμ = Broadcasting(Operation(⋅))(E,μ)
    φn = Broadcasting(Operation(×))(φ,n)
    Broadcasting(Operation(⋅))(φn,Eμ)
  end
  function gmom(φ,μ,ds) # Face gradient moment function: σ_F(φ,μ) = ∫((φ×n)·E∇μ)dF
    fmom(φ,Broadcasting(∇)(μ),ds)
  end
  function emom(φ,μ,ds) # Edge moment function: σ_E(φ,μ) = ∫((φ·t)*μ)dE
    t = get_edge_tangent(ds)
    φt = Broadcasting(Operation(⋅))(φ,t)
    Broadcasting(Operation(*))(φt,μ)
  end

  eb = LegendreBasis(Val(1), T, k, Polynomials._p_filter)  # 𝓟ₖ(e)
  moments = Tuple[(get_dimrange(p, 1), emom, eb)]
  for d in 2:D
    s = k + 1 - 2*(d-1)
    faces = get_dimrange(p, d)
    mom = d < D ? fmom : cmom
    if s > 0
      pb = LegendreBasis(Val(d), VectorValue{d,T}, s-1, Polynomials._p_filter) # 𝓟ₛ₋₁Λᵈ⁻¹
      push!(moments, (faces, mom, pb))
    end
    s < 0 && continue
    if d < D
      # the gradients of the Legendre polynomials of degree s+1 span ∇𝓗ₛ₊₁(f)
      # modulo 𝓟ₛ₋₁(f;R²), and are only permuted and sign-flipped by the
      # symmetries of the face
      gb = LegendreBasis(Val(d), T, s+1, Polynomials._p_complement_filter(s))
      push!(moments, (faces, gmom, gb))
    else
      gb = Polynomials._exact_poly_basis(Legendre, T, D, s, D-1; rotate_90=(D==2))
      push!(moments, (faces, cmom, gb))                                       # d𝓗ₛ₊₁Λᴰ⁻²
    end
  end
  conf = CurlConformity()
  dofs = MomentBasedDofBasis(p, prebasis, moments)
  face_own_dofs_permutations = moment_face_own_dofs_permutations(p, moments, conf)
  GenericRefFE{TrimmedSerendipityCurl}(
    length(dofs), p, prebasis, dofs, conf, nothing, get_face_own_moments(dofs);
    face_own_dofs_permutations)
end

function ReferenceFE(p::Polytope, ::TrimmedSerendipityCurl, ::Type{T}, order) where T
  TrimmedSerendipityCurlRefFE(T, p, order)
end

################################################################################
# Gluing on arbitrarily numbered meshes: the face DoFs of a cell are a signed
# permutation of the global ones, derived in closed form by
# `moment_face_own_dofs_permutations` (see MomentBasedReferenceFEs.jl) and
# stored in the reference FE; the signs are applied per cell by `FaceDofSignMap`
# (src/FESpaces/Pullbacks.jl), times the outward-normal signs where Gridap
# applies them. No basis of 𝓟ᵣ(□) is permuted without signs (the half turn has
# trace -1 on 𝓟₁), so signs are unavoidable.
