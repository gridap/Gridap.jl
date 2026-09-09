# RotatingPLambdaRefFEs.jl
#
# Reference FEs for the rotating P_rΛ¹ and trimmed P_r⁻Λ¹ bases: thin Gridap
# ReferenceFEs whose cross-cell H(curl) conformity is obtained through a
# per-cell change of basis given by the closed-form vertex-relabelling
# calculus of Gridap.Polynomials (RotatingPLambda/PΛRotations.jl and
# RotatingPLambda/PΛTrimmedRotations.jl). The per-cell change of basis itself
# (the compute_cell_bases_changes methods) lives in
# src/FESpaces/Pullbacks.jl.
#
# ─────────────────────────────────────────────────────────────────────────────
# DESIGN: virtually sorted cells
# ─────────────────────────────────────────────────────────────────────────────
# Every cell evaluates the SAME reference basis b (BarycentricPΛBasis or
# BarycentricPmΛBasis). H(curl) gluing across cells with ARBITRARY local vertex
# orderings is obtained by re-expressing each cell's covariant-Piola-pushed
# basis in the frame of the "virtually sorted" cell — the cell whose local
# vertices are relisted in increasing global-id order — via the rotation
# change of basis at
#
#   π_K = sortperm(cell global vertex ids)
#
# (sorted position → local index: π_K[i] is the local index holding the i-th
# smallest global vertex id). Both cells sharing a face then see the face in
# sorted global vertex order, so the face traces match DOF-by-DOF in bubble
# order: NO per-face data, NO sign flips, and identity face-own-dof
# permutations for every pindex.
#
# ─────────────────────────────────────────────────────────────────────────────
# DOF BASIS
# ─────────────────────────────────────────────────────────────────────────────
# The predofs are pointwise moments at BB lattice nodes, grouped per face:
# for a bubble entry (F, α, J) of face F,
#   • untrimmed, |α| = r: node x_β with β = α, covector Ψ_w;
#   • trimmed, |α| = r−1: node x_β with β = α + e_{J₁} (so the Whitney
#     covector is nonzero there), covector = the Whitney form ϕ^J evaluated
#     at x_β.
# The actual dof basis is the DUAL basis of the shape functions (= the basis
# itself), computed once on the reference element by the GenericRefFE
# predofs-constructor (dofs = compute_dofs(predofs, shapefuns)). This keeps
# shapefuns == prebasis == b, so the rotation calculus applies verbatim.


# ─────────────────────────────────────────────────────────────────────────────
# Names and singletons
# ─────────────────────────────────────────────────────────────────────────────

"""
    struct RotatingPΛName <: ReferenceFEName end

Reference FE name for the H(curl)-conforming vector proxied P_rΛ¹ element, supporting non-oriented meshes.
See [`RotatingPΛRefFE`](@ref). Its singleton instance is `rotating_pλ`.
"""
struct RotatingPΛName <: ReferenceFEName end

"""
    struct TrimmedPΛName <: ReferenceFEName end

Reference FE name for the H(curl)-conforming vector proxied P_r⁻Λ¹ element, supporting non-oriented meshes.
see [`TrimmedPΛRefFE`](@ref). Its singleton instance is `trimmed_pλ`.
"""
struct TrimmedPΛName  <: ReferenceFEName end

"""
    const rotating_pλ = RotatingPΛName()

Singleton of the [`RotatingPΛName`](@ref) reference FE name.
"""
const rotating_pλ = RotatingPΛName()

"""
    const trimmed_pλ = TrimmedPΛName()

Singleton of the [`TrimmedPΛName`](@ref) reference FE name.
"""
const trimmed_pλ  = TrimmedPΛName()

_pλ_basis(::RotatingPΛName, ::Type{T}, ::Val{D}, r, vertices; flavor) where {T,D} =
  BarycentricPΛBasis(Val(D), T, r, 1, vertices; flavor)
_pλ_basis(::TrimmedPΛName,  ::Type{T}, ::Val{D}, r, vertices; flavor) where {T,D} =
  BarycentricPmΛBasis(Val(D), T, r, 1, vertices; flavor)

Pushforward(::Type{RotatingPΛName}, ::CurlConformity) = CoVariantPiolaMap()
Pushforward(::Type{TrimmedPΛName},  ::CurlConformity) = CoVariantPiolaMap()

# ─────────────────────────────────────────────────────────────────────────────
# Predofs: pointwise moments at BB lattice nodes, grouped per face
# ─────────────────────────────────────────────────────────────────────────────
const _PΛBases = Union{BarycentricPΛBasis,BarycentricPmΛBasis}

# One pointwise moment per basis function, at the lattice node carrying it,
# contracted with the entry's own physical form direction.
function _pλ_predofs(b::_PΛBases, polytope, ::Val{D}) where D
  r      = get_order(b)
  V      = value_type(b)
  n_gf   = num_faces(polytope)
  verts  = get_vertex_coordinates(polytope)
  # bubble tables are keyed by sorted vertex sets, hence the sort
  fverts = map(v -> sort(collect(Int, v)), get_face_vertices(polytope))
  face_to_bfs = Dict(F => bfs for (F, bfs) in get_bubbles(b))

  # Exponent β (|β| = r) of the lattice node x_β = Σᵢ βᵢ/r · Vᵢ carrying the
  # predof of a bubble entry. The trimmed α has |α| = r−1, and is shifted onto
  # the first vertex of J so that the Whitney covector does not vanish at x_β.
  node_multiindex(::BarycentricPΛBasis, α, J) = copy(α)
  function node_multiindex(::BarycentricPmΛBasis, α, J)
    β = copy(α); β[J[1]] += 1; return β
  end

  # Predof covector at x_β: the entry's own (physical, Cartesian) form
  # direction evaluated at the node.  λ_i(x_β) = β_i/r for the Whitney form
  # ϕ^J = λ^{J₁} dλ^{J₂} − λ^{J₂} dλ^{J₁}, whose dλ^{J∖J(l)} are `b.m`.
  node_covector(b::BarycentricPΛBasis, w, β, J, sub_J_ids) = b.Ψ[w]
  node_covector(b::BarycentricPmΛBasis, w, β, J, sub_J_ids) =
    (β[J[1]]/r) * b.m[sub_J_ids[1]] - (β[J[2]]/r) * b.m[sub_J_ids[2]]

  lattice_point(β) = sum(β[i] * verts[i] for i in 1:D+1) / r

  all_nodes  = Point{D,Float64}[]
  f_moments  = Vector{Matrix{V}}(undef, n_gf)
  f_nodes    = Vector{UnitRange{Int}}(undef, n_gf)
  f_own_moms = Vector{Vector{Int}}(undef, n_gf)

  n_nodes, n_dofs  = 0, 0
  βs    = Vector{Vector{Int}}()
  β_row = Dict{Vector{Int},Int}()
  rows  = Int[]
  for gf in 1:n_gf
    bfs = get(face_to_bfs, fverts[gf], nothing)
    if isnothing(bfs) # face owns no bubble
      f_moments[gf]  = zeros(V, 0, 0)
      f_nodes[gf]    = 1:0
      f_own_moms[gf] = Int[]
      continue
    end
    n_bubble = length(bfs)
    resize!(βs, n_bubble)
    resize!(rows, n_bubble)
    empty!(β_row)
    for (j, (w, α, _, J)) in enumerate(bfs)
      β = node_multiindex(b, α, J)
      rows[j] = get!(β_row, β) do
        push!(all_nodes, lattice_point(β))
        n = length(β_row) + 1; βs[n] = β; return n
      end
    end
    n_βs = length(β_row)
    moms = zeros(V, (n_βs, n_bubble))
    for (j, (w, α, _, J, sub_J_ids)) in enumerate(bfs)
      moms[rows[j], j] = node_covector(b, w, βs[rows[j]], J, sub_J_ids)
    end
    f_moments[gf]  = moms
    f_nodes[gf]    = (n_nodes+1) : (n_nodes+n_βs)
    # The moment of a bubble entry is the dof dual to that entry, so the face
    # owns the entries' own indices w — which follow the bubble enumeration,
    # not the polytope face order this loop runs in.
    f_own_moms[gf] = collect(first(bfs)[1] : last(bfs)[1])
    n_nodes += n_βs
    n_dofs  += n_bubble
  end
  @assert n_dofs == length(b)

  return MomentBasedDofBasis(all_nodes, f_moments, f_nodes, f_own_moms)
end

# ─────────────────────────────────────────────────────────────────────────────
# Reference FE construction and the identity face-own-dof permutations
# ─────────────────────────────────────────────────────────────────────────────

function get_face_own_dofs_permutations(
    reffe::GenericRefFE{<:Union{RotatingPΛName,TrimmedPΛName}}, conf::Conformity)
  polytope      = get_polytope(reffe)
  face_own_dofs = get_face_own_dofs(reffe, conf)
  vtx_perms     = get_face_vertex_permutations(polytope)
  [[collect(1:length(own)) for _ in vtx_perms[gf]]
   for (gf, own) in enumerate(face_own_dofs)]
end

function _pλ_reffe(name::ReferenceFEName, ::Type{T}, p::Polytope{D}, r::Integer;
                   flavor=:BMM) where {T,D}
  @assert D in (2, 3) "only D = 2, 3 supported, got D = $D"
  @assert is_simplex(p) "only defined on simplices, got $p"
  @assert r ≥ 1 "r must be ≥ 1, got $r"
  basis     = _pλ_basis(name, T, Val(D), r, get_vertex_coordinates(p); flavor)
  n_dofs    = length(basis)
  predofs   = _pλ_predofs(basis, p, Val(D))
  face_dofs = get_face_own_funs(basis, p, CurlConformity())
  GenericRefFE{typeof(name)}(
    n_dofs, p, predofs, CurlConformity(), nothing, face_dofs, basis
  )
end

"""
    RotatingPΛRefFE(::Type{T}, p::Polytope{D}, r::Integer; flavor=:BMM)

H(curl)-conforming reference FE for the full P_rΛ¹ basis on the simplex `p`
(D = 2, 3; r ≥ 1). `T` is the type of scalar components of the vector proxied
shape function values. `flavor` selects the direction forms of the prebasis,
see [`BarycentricPΛBasis`](@ref).
"""
RotatingPΛRefFE(::Type{T}, p::Polytope, r; kwargs...) where T =
  _pλ_reffe(rotating_pλ, T, p, r; kwargs...)

"""
    TrimmedPΛRefFE(::Type{T}, p::Polytope{D}, r::Integer; flavor=:BMM)

H(curl)-conforming reference FE for the trimmed P_r⁻Λ¹ basis on the simplex `p`
(D = 2, 3; r ≥ 1). `T` is the type of scalar components of the vector proxied
shape function values. `flavor` selects the direction forms of the prebasis and,
for `:AFW`, makes its scalar factors Bernstein polynomials rather than bare
barycentric monomials, see [`BarycentricPmΛBasis`](@ref).
"""
TrimmedPΛRefFE(::Type{T}, p::Polytope, r; kwargs...) where T =
  _pλ_reffe(trimmed_pλ, T, p, r; kwargs...)

ReferenceFE(p::Polytope, ::RotatingPΛName, ::Type{T}, r; kwargs...) where T =
  RotatingPΛRefFE(T, p, r; kwargs...)
ReferenceFE(p::Polytope, ::TrimmedPΛName,  ::Type{T}, r; kwargs...) where T =
  TrimmedPΛRefFE(T, p, r; kwargs...)

