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
287, (8) and (12)]. 

## Prebasis

In 2D, implementation follows the augmented element approach of [Kirby, SMAI-JCM 4 (2018) 197]
We take `P₃(K;R²)`, of dimension 20, with constraints

    ∫_K (div Φ) q dK = 0    ∀ q ∈ P₂(K) ∩ P₀(K)^⊥
    ∫ₑ (Φ⋅n) μ ds    = 0    ∀ μ ∈ P₃(e) ∩ P₁(e)^⊥,  ∀ e

which enforce 11 constraints on the 20 prebasis functions, yielding the 9 DoFs.

In 3D, we directly build the prebasis `P₁(T;R³) + curl(b P₁(T;R³))`, of dimension 24. 

## Moments

Per facet `f` with unit normal `n`, `{q}` a basis of `P₁(f)` and `{w}` the rigid
motions of `f`:

    ℓ^{f,q}(Φ) = ∫_f (Φ⋅n) q,      ℓ^{f,w}(Φ) = ∫_f Φ_t ⋅ w

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
  cb = LegendreBasis(Val(1), T, 3, Polynomials._p_complement_filter(1))  # P₃(e) ∩ P₁(e)^⊥
  qb = DubinerBasis(Val(2), T, 2, Polynomials._p_complement_filter(0))   # P₂(K) ∩ P₀(K)^⊥
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
  fb = MonomialBasis(Val(2), T, 1, Polynomials._p_filter)                    # {1, u, v}
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

################################################################################
# Change of basis
#
# The cell-local map. The mesh-level `compute_cell_bases_changes`, which reads
# the orientation data off the model and maps this over the cells, lives in
# src/FESpaces/Pullbacks.jl.

# 2D: block diagonal, one 3x3 block per edge, acting on (ℓⁿ⁰, ℓⁿ¹, ℓᵗ⁰). Normal
# moments are preserved by the contravariant Piola map, while tangential ones
# pick up a normal component:
#
#   Wᵏ = [s 0 0; 0 s 0; α 0 β],  s = sign(det J),
#   α = (J n̂)⋅(J t̂)/|det J|,     β = ‖J t̂‖²/|det J|,
#
# post-multiplied by Dᵏ = diag(σᵏ,1,σᵏ). This is
# [Aznaran, Farrell & Kirby, SMAI-JCM 8 (2022) 399, (5.21)-(5.25)]; the
# degree-1 normal moment contributes an identity row, the parity of μ₁ being what
# keeps ℓⁿ¹ from mixing into ℓⁿ⁰.
struct MTWChangeOfBasis <: Map
  tangents::Vector{VectorValue{2,Float64}}
  normals::Vector{VectorValue{2,Float64}}
  transposed_inverse::Bool
end

function MTWChangeOfBasis(p::Polytope{2}, transposed_inverse::Bool)
  ts, ns = _edge_frames(p)
  MTWChangeOfBasis(ts, ns, transposed_inverse)
end

function return_cache(k::MTWChangeOfBasis, Jt, σ)
  ndofs = 3*length(k.tangents)
  CachedArray(zeros(Float64, ndofs, ndofs))
end

function evaluate!(cache, k::MTWChangeOfBasis, Jt, σ)
  nedges = length(k.tangents)
  setsize!(cache, (3*nedges, 3*nedges))
  M = cache.array
  fill!(M, zero(eltype(M)))

  s = sign(det(Jt))
  adetJ = abs(det(Jt))

  for e in 1:nedges
    # Jt is the transposed Jacobian, so `v ⋅ Jt` is `J v`.
    Jt̂ = k.tangents[e] ⋅ Jt
    Jn̂ = k.normals[e] ⋅ Jt
    α = (Jn̂ ⋅ Jt̂) / adetJ
    β = (Jt̂ ⋅ Jt̂) / adetJ
    σe = σ[e]

    # Dᵏ = diag(σ,1,σ) scales the columns of the block below.
    o = 3*(e - 1)
    if k.transposed_inverse
      # (Wᵏ)ᵀDᵏ
      M[o+1, o+1] = s * σe
      M[o+2, o+2] = s
      M[o+1, o+3] = α * σe
      M[o+3, o+3] = β * σe
    else
      # (Wᵏ)⁻¹Dᵏ
      M[o+1, o+1] = s * σe
      M[o+2, o+2] = s
      M[o+3, o+1] = -s * α * σe / β
      M[o+3, o+3] = σe / β
    end
  end

  return M
end

# 3D: one 6x6 block per face, and nothing else.
#
# The block is closed form. Its DoF weights are *equivariant* -- every one is a
# fixed reference object pushed forward by the map, `n = e₁×e₂` with `eᵢ = J êᵢ`,
# and the tangential weights `n×e₁`, `n×e₂`, `n×(x-a)` -- so the pulled-back
# physical weights expand in the reference weights with constant coefficients:
#
#   Jᵀwⱼ = [detJ·I₃  0 ;  Ã  λI₃] ŵ,   λ = |n|²/|n̂|²,  a = (Jn̂)⋅(n×e₁)/|n̂|²
#
# (and `b` likewise with `e₂`), whence `W = that / detJ`. The normal moments come
# out preserved exactly -- the identity block -- and `W⁻¹` is closed form too, so
# no quadrature, no cached shape-function moments and no matrix inversion appear
# at run time.
#
# This is the 3D analogue of [AFK, (5.21)-(5.25)], which they give only in 2D;
# FIAT/FInAT define the 3D element but transform it only on triangles. Getting it
# needs the equivariant weight basis: with an intrinsic frame (normalised
# tangents, unit normal) the coefficients stop being constant and the block no
# longer closes without integrating the shape functions.
struct TWChangeOfBasis <: Map
  face_dofs::Vector{Vector{Int}}          # dofs owned by each face
  # per (face, permutation), the reference edge vectors of the face taken in the
  # mesh's own vertex order -- the only reference data the block needs
  ref_edges::Matrix{NTuple{2,VectorValue{3,Float64}}}
  # per (face, permutation), `inv(T)` or `transpose(T)` as needed, where `T`
  # re-expresses the reference weights of the mesh-sorted frame in those of the
  # polytope's own frame -- see `_tw_frame_change`
  pre::Matrix{Matrix{Float64}}
  ndofs::Int
  transposed_inverse::Bool
end

function TWChangeOfBasis(reffe::ReferenceFE, transposed_inverse::Bool)
  p = get_polytope(reffe)
  own = get_face_own_dofs(reffe)
  fdim = get_dimrange(p, 2)
  face_dofs = [own[f] for f in fdim]
  nfaces = length(face_dofs)

  fverts = get_faces(p, 2, 0)
  vcoords = get_vertex_coordinates(p)
  vperms = get_face_vertex_permutations(p, 2)   # the list Gridap's pindex indexes

  # the same P⁻₁Λ¹ (RT₀) weight basis the tangential moments are declared with
  μb = FEEC_poly_basis(Val(2), Float64, 1, 1, :P⁻; rotate_90=true)
  nouts = get_facet_normal(p)                   # the outward normals `σ` uses
  nperms = length(first(vperms))
  ref_edges = Matrix{NTuple{2,VectorValue{3,Float64}}}(undef, nfaces, nperms)
  pre = Matrix{Matrix{Float64}}(undef, nfaces, nperms)
  for f in 1:nfaces, (pid, perm) in enumerate(vperms[f])
    # `pindex` maps cell-face-vertex -> global-face-vertex, so its inverse lists
    # the cell's local vertices in the order the mesh stores the face
    lv = fverts[f][invperm(perm)]
    â, b̂, ĉ = vcoords[lv[1]], vcoords[lv[2]], vcoords[lv[3]]
    ref_edges[f, pid] = (b̂ - â, ĉ - â)

    loc = fverts[f]
    T = _tw_frame_change((vcoords[loc[1]], vcoords[loc[2]], vcoords[loc[3]]),
                         (â, b̂, ĉ), μb, perm, nouts[f])
    pre[f, pid] = transposed_inverse ? transpose(T) : inv(T)
  end

  TWChangeOfBasis(face_dofs, ref_edges, pre, num_dofs(reffe), transposed_inverse)
end

# The 6x6 matrix `T` with `ŵˢⱼ = Σₖ Tⱼₖ ŵˡₖ`, expressing the canonical weights of
# the mesh-sorted frame of a face in the ones the moments actually declare, in
# the polytope's own frame. Both blocks are closed form.
#
# The normal block is `ε A`, with `ε = ±1` the sign the scaled normal `n̂ₛ = ê₁×ê₂`
# picks up under the relabelling and `A` the affine change from `(1, uˢ, vˢ)` to
# `(1, uˡ, vˡ)`, determined exactly by three points.
#
# The tangential block is `Mμ⁻¹ R`, where `R = rotation_change_of_basis(μb,
# invperm(π))` is the vertex relabelling of `μ` and `Mμ` holds the `RT₀`
# coefficients `(a, b, c)` of `μⱼ = (a + c u, b + c v)`, which is the change from
# `{Ĵf μⱼ}` to the canonical `{ê₁, ê₂, x̂-â}`. `μ` is a trimmed FEEC space, so `R`
# is available in closed form and, at this degree, is a signed permutation with
# integer entries. See Badia, Manyer & Marteau, *Rotating bases for finite
# element exterior calculus*.
#
# Both factors are exact, so nothing here is fitted. `Mμ` is geometry-free and
# `R` is integer; the only reason `T` is assembled per (face, permutation) at all
# is the normal block's dependence on the face's own vertex labelling.
function _tw_frame_change(lverts, sverts, μb, π, nout)
  v3(w) = Float64[w[1], w[2], w[3]]
  T = zeros(6, 6)

  a, e1, e2 = lverts[1], lverts[2] - lverts[1], lverts[3] - lverts[1]
  as, es1, es2 = sverts[1], sverts[2] - sverts[1], sverts[3] - sverts[1]
  ns = cross(e1, e2)
  ε = ns ⋅ cross(es1, es2) > 0 ? 1.0 : -1.0
  # the moments state their weights with the polytope's *outward* normal, which
  # carries no orientation of the shared face. `sout` is the sign relating it to
  # ê₁×ê₂, the orientation the sorted frame does carry, and `area2` undoes the
  # extension's normalisation of the tangential weights. Both are per-face
  # constants.
  sout = nout ⋅ ns > 0 ? 1.0 : -1.0
  area2 = norm(ns)
  Es = hcat(v3(es1), v3(es2))

  # (1, uˢ, vˢ) against (1, uˡ, vˡ) at three points fixes the affine relabelling
  Vl = zeros(3, 3)
  Vs = zeros(3, 3)
  for (i, (u, v)) in enumerate(((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)))
    uv = Es \ v3((a + u * e1 + v * e2) - as)
    Vl[i, :] = [1.0, u, v]
    Vs[i, :] = [1.0, uv[1], uv[2]]
  end
  T[1:3, 1:3] = sout * ε * transpose(Vl \ Vs)

  # μⱼ = (a + c u, b + c v), so Ĵf μⱼ = a ê₁ + b ê₂ + c (x̂ - â)
  m00 = evaluate(μb, [Point(0.0, 0.0)])
  m10 = evaluate(μb, [Point(1.0, 0.0)])
  Mμ = zeros(3, 3)
  for j in 1:3
    Mμ[j, :] = [m00[1, j][1], m00[1, j][2], m10[1, j][1] - m00[1, j][1]]
  end
  T[4:6, 4:6] = sout * area2 *
                inv(Mμ) * Polynomials.rotation_change_of_basis(μb, invperm(collect(π)))
  T
end

function return_cache(k::TWChangeOfBasis, Jt, pids)
  CachedArray(zeros(Float64, k.ndofs, k.ndofs)), zeros(6, 6)
end

function evaluate!(_cache, k::TWChangeOfBasis, Jt, pids)
  cache, W = _cache
  setsize!(cache, (k.ndofs, k.ndofs))
  M = cache.array
  fill!(M, zero(eltype(M)))
  detJt = det(Jt)
  adetJ = abs(detJt)
  sgn = sign(detJt)

  for f in eachindex(k.face_dofs)
    ê1, ê2 = k.ref_edges[f, pids[f]]
    n̂ = cross(ê1, ê2)
    s2 = n̂ ⋅ n̂
    e1 = ê1 ⋅ Jt                        # `v ⋅ Jt` is `J v`
    e2 = ê2 ⋅ Jt
    n = cross(e1, e2)
    Jn̂ = n̂ ⋅ Jt

    μ = (n ⋅ n) / (s2 * adetJ)
    a = (Jn̂ ⋅ cross(n, e1)) / (s2 * adetJ)
    b = (Jn̂ ⋅ cross(n, e2)) / (s2 * adetJ)

    fill!(W, 0.0)
    if k.transposed_inverse                       # Wᵀ
      for i in 1:3
        W[i, i] = sgn
        W[3+i, 3+i] = μ
      end
      W[1, 4] = a; W[1, 5] = b; W[2, 6] = a; W[3, 6] = b
    else                                          # W⁻¹
      iμ = 1 / μ
      for i in 1:3
        W[i, i] = sgn
        W[3+i, 3+i] = iμ
      end
      W[4, 1] = -sgn * a * iμ; W[5, 1] = -sgn * b * iμ
      W[6, 2] = -sgn * a * iμ; W[6, 3] = -sgn * b * iμ
    end

    # the physical dofs are stated in the mesh-sorted frame, the reference ones
    # in the polytope's own, so the frame change composes in
    P = k.pre[f, pids[f]]
    dofs = k.face_dofs[f]
    for i in 1:6, j in 1:6
      M[dofs[i], dofs[j]] = sum(P[i, m] * W[m, j] for m in 1:6)
    end
  end

  return M
end
