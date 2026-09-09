############################################################################################
# The Mardal-Tai-Winther family, on triangles and on tetrahedra.
#
# One element in two dimensions:
#
#   V(K) = P₁(K;Rᴰ) + curl(b P₁(K;Λ^{D-2})),    b the barycentric bubble
#
# with `curl` the scalar rot in 2D and the vector curl in 3D. `b` is cubic in 2D
# and quartic in 3D, and `curl` drops one degree, so `V(K)` sits in `P₃` and `P₄`
# respectively, of dimension 9 and 24.
#
# The DoFs are the same in both dimensions, six per facet in 3D and three per
# facet in 2D: the normal component against `P₁(f)`, and the tangential component
# against the *rigid motions* of `f`. The rigid motions of a segment are the
# constants, one-dimensional, which is why 2D has a single tangential moment per
# edge; those of a triangle are three-dimensional, two translations and the
# in-plane rotation. FIAT unifies the two the same way.
#
# The element is H(div)-conforming and only weakly H¹ -- the tangential jump
# across a facet does not vanish, but its rigid-motion moments do -- which is
# what makes it robust for Darcy-Stokes uniformly in the parameter.
#
# THE TWO CONSTRUCTIONS ARE NOT THE SAME
#
# 2D builds `V(K)` as a *kernel*: [Mardal, Tai & Winther, §4.1] define it by
# `div Φ ∈ P₀(K)` and `(Φ⋅n)|ₑ ∈ P₁(e)`, eleven constraints inside the
# 20-dimensional `P₃(K;R²)`, and the element is assembled as the augmented
# element of [Kirby, SMAI-JCM 4 (2018) 197, §5.5].
#
# 3D builds it as a *sum*: `P₁ + curl(bP₁)` has no comparable constraint set, so
# the 24 generators are written into the Bernstein basis of `P₄(T;R³)` by exact
# integer arithmetic instead.
#
# The 2D space is in fact the same sum -- verified -- so the 2D case could be
# rebuilt the 3D way, which would drop the constraints, the `restrict` and the
# 20x20 inverse. Left as is for now: the 2D change of basis is the published
# [Aznaran, Farrell & Kirby, SMAI-JCM 8 (2022) 399, (5.21)-(5.25)] and is worth
# keeping in the form it was verified in.
#
# Consequently the change of basis has two shapes: 3x3 per edge in 2D, keyed by
# an edge orientation sign, and 6x6 per face in 3D, keyed by the face's `pindex`
# -- a triangular face has 3! = 6 orderings, so the mismatch between neighbours
# is a permutation and not a sign. Both are closed form. See `FESpaces`.
############################################################################################

"""
    struct MardalTaiWinther <: ReferenceFEName end

Reference FE name for the Mardal--Tai--Winther element, H(div)-conforming and
weakly H¹, on triangles and tetrahedra. See [`MardalTaiWintherRefFE`](@ref). Its
singleton instance is [`mtw`](@ref).
"""
struct MardalTaiWinther <: ReferenceFEName end

"""
    const mtw = MardalTaiWinther()

Singleton of the [`MardalTaiWinther`](@ref) reference FE name.
"""
const mtw = MardalTaiWinther()

Pushforward(::Type{MardalTaiWinther}) = ContraVariantPiolaMap()

"""
    MardalTaiWintherRefFE(::Type{T}, p::Polytope{D})

The Mardal--Tai--Winther reference FE on the simplex `K`, with `T` the scalar
type, for `D = 2` and `D = 3`:

    MTW(K) = P₁(K;Rᴰ) + curl(b P₁(K;Λ^{D-2})),    b the barycentric bubble

of dimension 9 on a triangle [Mardal, Tai & Winther, SIAM J. Numer. Anal. 40
(2002) 1605, Lemma 4.1] and 24 on a tetrahedron [Tai & Winther, Calcolo 43 (2006)
287, (8) and (12)]. `curl` is the scalar rot in 2D and the vector curl in 3D, so
the space is cubic on a triangle and quartic on a tetrahedron. It is
H(div)-conforming and weakly H¹: the tangential jump across a facet does not
vanish, but its moments against the rigid motions of that facet do.

## Prebasis

In 2D, `P₃(K;R²)` of dimension 20 — the *ambient* space, not `MTW(K)` itself,
which is cut out of it by the constraints below and assembled as the augmented
element of [Kirby, SMAI-JCM 4 (2018) 197, §5.5]; `length(get_prebasis(reffe))` is
therefore larger than `num_dofs(reffe)` there.

In 3D, `MTW(K)` itself, of dimension 24. The space is a *sum* and not the kernel
of a constraint set, so the augmented construction does not apply and the 24
generators are written into the Bernstein basis of `P₄(K;R³)` directly.

## Moments

Per facet `f` with unit normal `n`, `{q}` a basis of `P₁(f)` and `{w}` the rigid
motions of `f`:

    ℓ^{f,q}(Φ) = ∫_f (Φ⋅n) q,      ℓ^{f,w}(Φ) = ∫_f Φ_t ⋅ w

The rigid motions of an edge are the constants, so 2D has two normal moments and
one tangential moment per edge, `3·3 = 9`; those of a triangle are the two
constant tangent fields and the in-plane rotation, so 3D has three of each per
face, `6·4 = 24`. The first group is what gives H(div) conformity, the second the
weak continuity.

In 3D the rigid motions are those of the *physical* face, taken in each face's
own vertices in increasing global id order. An affine extension of the reference
triangle's rigid motions does not carry its rotation to the face's own, and with
those weights the reference and physical DoFs span different spaces of
functionals — measurably so: the generalized Vandermonde stops being block
diagonal.

## Constraints

In 2D only, the eleven that cut `MTW(K)` out of `P₃(K;R²)`:

    ∫_K (div Φ) q dK = 0    ∀ q ∈ P₂(K) ∩ P₀(K)^⊥            (5)
    ∫ₑ (Φ⋅n) μ ds    = 0    ∀ μ ∈ P₃(e) ∩ P₁(e)^⊥,  ∀ e      (2 per edge)

`div Φ` lies in `P₂(K)` a priori and is pinned to its `P₀` part, and `(Φ⋅n)|ₑ`
lies in `P₃(e)` a priori and is pinned to its `P₁` part. `9 + 11 = 20`.
"""
function MardalTaiWintherRefFE(::Type{T}, p::Polytope{D}) where {T,D}
  # `@notimplemented`, not `@check`: the latter is compiled out in Gridap's
  # performance mode, and this guard must fire in every mode.
  @notimplementedif !(D in (2, 3) && is_simplex(p)) """\n
  The Mardal-Tai-Winther element is only defined on 2D and 3D simplices, got $p.
  """
  _mtw_reffe(T, p)
end

function ReferenceFE(p::Polytope, ::MardalTaiWinther, ::Type{T}) where T
  MardalTaiWintherRefFE(T, p)
end

function ReferenceFE(p::Polytope, ::MardalTaiWinther, ::Type{T}, order) where T
  @notimplementedif order != 1 """\n
  The Mardal-Tai-Winther element exists for order 1 only, got $order.
  """
  MardalTaiWintherRefFE(T, p)
end

# Identity for every admissible vertex permutation of every face. In 2D reversing
# an edge changes the sign of two of its DoFs but never their order; in 3D the
# physical DoFs are defined from each face's vertices in increasing global id
# order, so both neighbours agree on them, in the same order, whatever order they
# list the face themselves. Either way the change of basis carries the rest.
function get_face_own_dofs_permutations(
  reffe::GenericRefFE{MardalTaiWinther}, conf::Conformity
)
  _identity_dof_permutations(reffe, conf)
end

# ─────────────────────────────────────────────────────────────────────────────
# 2D: the augmented element
#
# The 9 DoFs and the 11 constraints are *both* declared as ordinary moments on
# the unconstrained Bernstein basis of P₃(K;R²). Together they are a dual basis
# for all of it -- a Φ killed by all twenty is in MTW(K) by the constraints and
# then zero by unisolvency -- so inverting their 20x20 Vandermonde and keeping
# the columns dual to the DoFs gives the nodal basis.
#
# The constraints need two `MomentBasedDofBasis`es because a single one carries
# one `operator` for all its moments, and the divergence rows need `∇` while
# everything else needs none. `vcat` joins them.
#
# The weight bases are Legendre on the edges and Dubiner on the cell. Both are
# L²-orthogonal on the face they live on, so both `P₃(e) ∩ P₁(e)^⊥` and
# `P₂(K) ∩ P₀(K)^⊥` are a plain selection by degree -- `_p_complement_filter`,
# no projection. Legendre cannot serve on the cell: it is a tensor product, so
# orthogonal on a segment and an n-cube but not on a triangle.
# ─────────────────────────────────────────────────────────────────────────────

function _mtw_reffe(::Type{T}, p::Polytope{2}) where T
  prebasis = BernsteinBasisOnSimplex(Val(2), VectorValue{2,T}, 3)

  # DoF moments
  nb = LegendreBasis(Val(1), T, 1)   # μ₀, μ₁ : the normal DoF weights
  tb = LegendreBasis(Val(1), T, 0)   # μ₀     : the tangential DoF weight
  function nmom(φ, μ, ds)            # ∫ₑ (Φ⋅n) μ ds
    n = _edge_normal(ds)
    Broadcasting(Operation(*))(Broadcasting(Operation(⋅))(φ, n), μ)
  end
  function tmom(φ, μ, ds)            # ∫ₑ (Φ⋅t) μ ds
    t = get_edge_tangent(ds)
    Broadcasting(Operation(*))(Broadcasting(Operation(⋅))(φ, t), μ)
  end

  # Constraint moments. The normal-trace constraints share `nmom` and so ride on
  # the same (identity-operator) basis as the DoFs; the divergence ones need `∇`.
  cb = LegendreBasis(Val(1), T, 3, _p_complement_filter(1))  # P₃(e) ∩ P₁(e)^⊥
  qb = DubinerBasis(Val(2), T, 2, _p_complement_filter(0))   # P₂(K) ∩ P₀(K)^⊥
  Id = ConstantField(TensorValue(one(T), zero(T), zero(T), one(T)))
  divmom(φ, μ, ds) = Broadcasting(Operation(*))(   # ∫_K (∇Φ ⊙ I) μ dK = ∫_K (div Φ) μ dK
    Broadcasting(Operation(⊙))(φ, Id), μ
  )

  edges = get_dimrange(p, 1)
  cell = get_dimrange(p, 2)
  moments = Tuple[
    (edges, nmom, nb), (edges, tmom, tb),                      # Edge DoFs
    (edges, nmom, cb),                                         # Edge constraints
  ]
  constraints = Tuple[(cell, divmom, qb)]                      # Cell constraints

  edge_functionals = MomentBasedDofBasis(p, prebasis, moments)
  div_functionals = MomentBasedDofBasis(p, prebasis, constraints, ∇)
  full = vcat(edge_functionals, div_functionals)

  # each edge contributes (ℓⁿ⁰, ℓⁿ¹, ℓᵗ⁰, c₂, c₃) and the first three are DoFs
  nedges = num_faces(p, 1)
  dof_ids = reduce(vcat, [5*(k - 1) .+ (1:3) for k in 1:nedges])
  ndofs = length(dof_ids)

  @check length(full) == length(prebasis) """\n
  The augmented element needs as many functionals as prebasis functions, got
  $(length(full)) and $(length(prebasis)).
  """
  shapefuns = linear_combination(inv(evaluate(full, prebasis))[:, dof_ids], prebasis)
  dofs = restrict(edge_functionals, dof_ids)

  # the post-condition of the whole construction, and the thing that catches a
  # rank-deficient or badly conditioned constraint set
  @check maximum(abs, evaluate(dofs, shapefuns) - Matrix{Float64}(I, ndofs, ndofs)) < 1e-10 """\n
  The augmented MTW Vandermonde is singular or badly conditioned.
  """
  GenericRefFE{MardalTaiWinther}(
    ndofs, p, prebasis, dofs, DivConformity(), nothing,
    get_face_own_moments(dofs), shapefuns
  )
end

# ─────────────────────────────────────────────────────────────────────────────
# 3D: the explicit sum
# ─────────────────────────────────────────────────────────────────────────────

# `V(T) = P₁(T;R³) + curl(b P₁(T;R³))`, expressed in the Bernstein basis of
# `P₄(p;R³)` -- which contains it, `b·λ` being quintic and its curl quartic.
#
# In barycentric coordinates this is exact arithmetic and two textbook
# identities, because the element is *defined* barycentrically:
#
#   * `P₁(T;R³) = span{λᵢ eⱼ}`, and degree elevation gives
#     `λᵢ = Σ_α (αᵢ/K) B_α⁽ᴷ⁾`;
#   * `b·λ_m = λ₀λ₁λ₂λ₃λ_m` is, up to a constant, the *single* Bernstein
#     function `B_γ⁽⁵⁾` with `γ = (1,1,1,1) + e_m` -- no product to expand;
#   * `curl(f e) = ∇f × e` for constant `e`, and
#     `∇B_γ⁽ᴷ⁾ = K Σᵢ B_{γ-eᵢ}⁽ᴷ⁻¹⁾ ∇λᵢ` lands the result directly in the
#     degree-4 basis, with `∇λᵢ × eⱼ` constant.
#
# Every `γᵢ ≥ 1`, so no term of the gradient sum drops and nothing needs
# guarding. Only `rank(C) == 24` is asserted: containment in `P₄` holds by
# construction. Bernstein rather than monomials because this basis is *why* the
# construction is clean, and because it is far better conditioned on a simplex.
function _mtw_prebasis_3d(::Type{T}, p::Polytope{3}) where T
  @notimplementedif get_vertex_coordinates(p) != [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
                                                  Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0)] """\n
  The 3D Mardal-Tai-Winther prebasis assumes Gridap's reference tetrahedron,
  whose barycentric gradients ∇λ are hard-coded below.
  """
  amb = BernsteinBasisOnSimplex(Val(3), VectorValue{3,T}, 4)
  terms = bernstein_terms(4, 3)
  id = Dict(Tuple(a) => k for (k, a) in enumerate(terms))
  # a MultiValue basis scatters the components of each term consecutively
  row(α, c) = 3*(id[α] - 1) + c

  ∇λ = (VectorValue{3,T}(-1, -1, -1), VectorValue{3,T}(1, 0, 0),
        VectorValue{3,T}(0, 1, 0),    VectorValue{3,T}(0, 0, 1))
  es = (VectorValue{3,T}(1, 0, 0), VectorValue{3,T}(0, 1, 0), VectorValue{3,T}(0, 0, 1))

  C = zeros(T, length(amb), 24)
  k = 0
  for i in 1:4, j in 1:3                      # λᵢ eⱼ, raised to degree 4
    k += 1
    for (a, α) in enumerate(terms)
      C[3*(a - 1) + j, k] = α[i] / 4
    end
  end
  for m in 1:4, j in 1:3                      # curl(B_γ⁽⁵⁾ eⱼ) = ∇B_γ⁽⁵⁾ × eⱼ
    k += 1
    γ = ntuple(q -> q == m ? 2 : 1, 4)
    for i in 1:4
      v = cross(∇λ[i], es[j])
      β = ntuple(q -> q == i ? γ[q] - 1 : γ[q], 4)
      for c in 1:3
        C[row(β, c), k] += 5 * v[c]
      end
    end
  end

  @check rank(C) == 24 """\n
  The 3D Mardal-Tai-Winther space should be 24-dimensional, got rank $(rank(C)).
  """
  linear_combination(C, amb)
end

# The 24 DoFs, six per face: three moments of the normal component against
# `P₁(f)`, then three of the tangential component against the rigid motions of
# the face. Declared through the `(faces, σ, μ)` moment API with the stock
# accessors: `get_facet_normal(ds)` for the outward normal and
# `get_extension(ds)` for the facet's contravariant Piola map, which *preserves*
# `RT₀` and so lets the tangential weights be stated through `μ`. The outward
# normal carries no orientation of the shared face; that is supplied in the
# change of basis.
#
# `fb` is `{1, u, v}` on the reference triangle, the weight of the normal
# moments. `fb_rt0` is `RT₀ = P⁻₁Λ¹`, three-dimensional, which is what the rigid
# motions are -- `RM(f)` is 3-dimensional at any degree and coincides with RT
# only at lowest order, so it is not a family parameter. FIAT makes the same
# choice.
function _mtw_dof_basis_3d(::Type{T}, p::Polytope{3}, prebasis) where T
  fb = MonomialBasis(Val(2), T, 1, _p_filter)                    # {1, u, v}
  fb_rt0 = FEEC_poly_basis(Val(2), T, 1, 1, :P⁻; rotate_90=true) # P⁻₁Λ¹ on the face

  function nmom(φ, μ, ds)                     # ∫_f (v⋅n) q dA,  q ∈ P₁(f)
    n = get_facet_normal(ds)
    Broadcasting(Operation(*))(Broadcasting(Operation(⋅))(φ, n), μ)
  end
  function rmom(φ, μ, ds)                     # ∫_f v ⋅ (n × Ĵf μ) dA,  μ ∈ RT₀(f̂)
    n = get_facet_normal(ds)
    w = Broadcasting(Operation(⋅))(get_extension(ds), μ)
    Broadcasting(Operation(⋅))(φ, Broadcasting(Operation(cross))(n, w))
  end

  faces = get_dimrange(p, 2)
  MomentBasedDofBasis(p, prebasis, Tuple[(faces, nmom, fb), (faces, rmom, fb_rt0)])
end

function _mtw_reffe(::Type{T}, p::Polytope{3}) where T
  prebasis = _mtw_prebasis_3d(T, p)
  dofs = _mtw_dof_basis_3d(T, p, prebasis)
  face_own_dofs = get_face_own_moments(dofs)
  @check length(dofs) == length(prebasis) == 24

  GenericRefFE{MardalTaiWinther}(
    24, p, prebasis, dofs, DivConformity(), nothing, face_own_dofs
  )
end
