"""
    struct Argyris <: ReferenceFEName end

Reference FE name for the quintic Argyris triangle, a C¹-conforming element. See
[`ArgyrisRefFE`](@ref). Its singleton instance is [`argyris`](@ref).
"""
struct Argyris <: ReferenceFEName end

"""
    const argyris = Argyris()

Singleton of the [`Argyris`](@ref) reference FE name.
"""
const argyris = Argyris()

# Mapped by the plain pullback; this is also the default, stated for clarity.
Pushforward(::Type{Argyris}) = IdentityPiolaMap()

function _argyris_dof_basis(::Type{T}, p::Polytope{2}, prebasis) where T
  c0 = MonomialBasis(Val(0), T, 0)  # the constant on a vertex
  e0 = MonomialBasis(Val(1), T, 0)  # the constant on an edge
  verts = get_dimrange(p, 0)
  edges = get_dimrange(p, 1)

  # Directions contracted with ∇u and ∇∇u to pick single components. `ConstantField`
  # is what lifts a fixed value into a `Field`, which `Operation` needs; a bare
  # `VectorValue` or `TensorValue` has no `evaluate!`.
  Ei = (ConstantField(VectorValue(one(T), zero(T))),      # e₁
        ConstantField(VectorValue(zero(T), one(T))))      # e₂
  Eij = (ConstantField(TensorValue(one(T), zero(T), zero(T), zero(T))),   # e₁⊗e₁
         ConstantField(TensorValue(zero(T), zero(T), one(T), zero(T))),   # e₁⊗e₂
         ConstantField(TensorValue(zero(T), zero(T), zero(T), one(T))))   # e₂⊗e₂
  Emom(E) = (φ, μ, ds) -> Broadcasting(Operation(*))(     # ∇u(v)⊙E or ∇∇u(v)⊙E
    Broadcasting(Operation(⊙))(φ, E), μ
  )

  vmom(φ, μ, ds) = Broadcasting(Operation(*))(φ, μ)
  function emom(φ, μ, ds)  # σ_e(∇u,μ) = ∫ₑ (∇u⋅n) μ ds
    n = _edge_normal(ds)
    Broadcasting(Operation(*))(Broadcasting(Operation(⋅))(φ, n), μ)
  end

  values = MomentBasedDofBasis(p, prebasis, Tuple[(verts, vmom, c0)])
  grads = MomentBasedDofBasis(
    p, prebasis, Tuple[(verts, Emom(E), c0) for E in Ei], ∇
  )
  hessians = MomentBasedDofBasis(
    p, prebasis, Tuple[(verts, Emom(E), c0) for E in Eij], ∇∇
  )
  normals = MomentBasedDofBasis(p, prebasis, Tuple[(edges, emom, e0)], ∇)

  return vcat(values, grads, hessians, normals)
end

function _argyris_face_own_dofs(p::Polytope{2})
  nv, ne = num_faces(p, 0), num_faces(p, 1)
  own = [Int[] for _ in 1:num_faces(p)]
  for v in 1:nv
    own[v] = [v, nv + 2*v - 1, nv + 2*v, 3*nv + 3*v - 2, 3*nv + 3*v - 1, 3*nv + 3*v]
  end
  for e in 1:ne
    own[nv+e] = [6*nv + e]
  end
  return own
end

"""
    ArgyrisRefFE(::Type{T}, p::Polytope{2})

The Argyris reference FE on the triangle `K`, with `T` the scalar type: the
quintic C¹ plate element of [Argyris, Fried & Scharpf, Aero. J. 72 (1968) 701],
21 DoFs, with

    Argyris(K) = P₅(K)

of dimension 21. Unconstrained; what is non-standard is the DoF set, which
contains point-evaluated *derivatives* and so is not preserved by the plain
pullback — the element therefore needs a cell-dependent change of basis, see
`compute_cell_bases_changes` in `FESpaces`.

## Prebasis

`P₅(K)`, of dimension 21, the element's space itself.

## Moments

At each vertex `v` the value, the two gradient components and the three
independent Hessian components; per edge `e` with unit tangent `t` and normal
`n = R t`, the mean normal derivative:

    ℓ^{v,0}(u)   = u(v)
    ℓ^{v,i}(u)   = ∇u(v) ⋅ eᵢ,               i = 1, 2
    ℓ^{v,ij}(u)  = ∇∇u(v) ⊙ E_ij,            E_ij ∈ {e₁⊗e₁, e₁⊗e₂, e₂⊗e₂}
    ℓ^e(u)       = ∫ₑ (∇u⋅n) ds

`6·3 + 3 = 21`. As in [`MorleyRefFE`](@ref) the edge DoF is a moment, not the
classical `∂u/∂n(m_e)`; the two give the *same* global C¹ space, since on an edge
`∂u/∂n` is a quartic four of whose coefficients are already fixed by the shared
vertex DoFs, so either choice pins the fifth. The moment is preferred because it
makes the node completion of [Kirby, SMAI-JCM 4 (2018) 197, §3.3.2] collapse to
`u(v_b) - u(v_a)`, exact for any smooth `u`, where the midpoint form needs six
terms per edge and is exact only on `P₅`.

The `E_ij` are written out rather than taken from the component basis of the
symmetric matrices, so that the DoFs are exactly `∂ᵢu` and `∂ᵢⱼu`: contracting
two symmetric matrices sums over both off-diagonal slots, which would make the
mixed second derivative `2∂₁₂u`.

The conformity is `H1Conformity()`, which in Gridap fixes *which faces own and
share DoFs* — here all of them. The assembled space is in fact C¹; Gridap has no
conformity object saying so, and none is needed, since the gluing pattern is what
the machinery consumes.
"""
function ArgyrisRefFE(::Type{T}, p::Polytope{D}) where {T,D}
  @notimplementedif !(D == 2 && is_simplex(p)) """\n
  The Argyris element is only defined on 2D simplices, got $p.
  """
  prebasis = BernsteinBasisOnSimplex(Val(2), T, 5)  # P₅(K), dim 21
  dofs = _argyris_dof_basis(T, p, prebasis)
  face_own_dofs = _argyris_face_own_dofs(p)

  ndofs = 6*num_faces(p, 0) + num_faces(p, 1)
  @check length(dofs) == ndofs == length(prebasis)

  GenericRefFE{Argyris}(ndofs, p, prebasis, dofs, H1Conformity(), nothing, face_own_dofs)
end

function ReferenceFE(p::Polytope, ::Argyris, ::Type{T}) where T
  ArgyrisRefFE(T, p)
end

function ReferenceFE(p::Polytope, ::Argyris, ::Type{T}, order) where T
  @notimplementedif order != 5 """\n
  The Argyris element exists for order 5 only, got $order.
  """
  ArgyrisRefFE(T, p)
end

function get_face_own_dofs_permutations(reffe::GenericRefFE{Argyris}, conf::Conformity)
  _identity_dof_permutations(reffe, conf)
end
