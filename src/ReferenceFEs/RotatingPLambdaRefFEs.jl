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
# Every cell evaluates the SAME reference basis b (RotatingPΛBasis or
# TrimmedPΛBasis). H(curl) gluing across cells with ARBITRARY local vertex
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
# for a bubble entry of face F,
#   • untrimmed (F, k, α), |α| = r: node x_β with β = α, covector ψ(F,k,α);
#   • trimmed (F, (e1,e2), α), |α| = r−1: node x_β with β = α + e_{e1} (so
#     the Whitney covector is nonzero there), covector = the Whitney form
#     ϕ(e1,e2) evaluated at x_β.
# The actual dof basis is the DUAL basis of the shape functions (= the basis
# itself), computed once on the reference element by the GenericRefFE
# predofs-constructor (dofs = compute_dofs(predofs, shapefuns)). This keeps
# shapefuns == prebasis == b, so the rotation calculus applies verbatim.


# ─────────────────────────────────────────────────────────────────────────────
# Names and singletons
# ─────────────────────────────────────────────────────────────────────────────

"""
    struct RotatingPΛName <: ReferenceFEName end

Reference FE name for the H(curl)-conforming rotating P_rΛ¹ element, see
[`RotatingPΛRefFE`](@ref). Its singleton instance is `rotating_pλ`.
"""
struct RotatingPΛName <: ReferenceFEName end

"""
    struct TrimmedPΛName <: ReferenceFEName end

Reference FE name for the H(curl)-conforming trimmed P_r⁻Λ¹ element, see
[`TrimmedPΛRefFE`](@ref). Its singleton instance is `trimmed_pλ`.
"""
struct TrimmedPΛName  <: ReferenceFEName end

const rotating_pλ = RotatingPΛName()
const trimmed_pλ  = TrimmedPΛName()

_pλ_basis(::RotatingPΛName, ::Val{D}, r) where D = RotatingPΛBasis(Val(D), Float64, r)
_pλ_basis(::TrimmedPΛName,  ::Val{D}, r) where D = TrimmedPΛBasis(Val(D), Float64, r)

# Both bases are H(curl) conforming 1-form spaces: covariant Piola map.
Pushforward(::Type{RotatingPΛName}, ::CurlConformity) = CoVariantPiolaMap()
Pushforward(::Type{TrimmedPΛName},  ::CurlConformity) = CoVariantPiolaMap()

# ─────────────────────────────────────────────────────────────────────────────
# Covariant Piola map on DifferentialFormValue 1-forms
# ─────────────────────────────────────────────────────────────────────────────

# Pushforward: ω_phys = J⁻ᵀ ω_ref
@inline function evaluate!(cache, ::CoVariantPiolaMap,
  v::DifferentialFormValue{1,D,T,D,Cartesian{D}}, Jt::Number) where {D,T}
  vv = VectorValue{D,T}(v.data)
  DifferentialFormValue{1,D}(Tuple((pinvJt(Jt) ⋅ vv).data))
end

# Inverse pushforward: ω_ref = Jᵀ ω_phys
@inline function evaluate!(cache, ::InversePushforward{CoVariantPiolaMap},
  v::DifferentialFormValue{1,D,T,D,Cartesian{D}}, Jt::Number) where {D,T}
  vv = VectorValue{D,T}(v.data)
  DifferentialFormValue{1,D}(Tuple((Jt ⋅ vv).data))
end

# ─────────────────────────────────────────────────────────────────────────────
# Predofs: pointwise moments at BB lattice nodes, grouped per face
# ─────────────────────────────────────────────────────────────────────────────

# One pointwise moment per basis function, at the lattice node carrying it,
# contracted with the entry's own physical form direction.
function _pλ_predofs(b::_PΛBases, polytope, ::Val{D}) where D
  r      = get_order(b)
  V      = value_type(b)
  n_gf   = num_faces(polytope)
  verts  = get_vertex_coordinates(polytope)
  # bubble tables are keyed by sorted vertex sets, hence the sort
  fverts = map(v -> sort(collect(Int, v)), get_face_vertices(polytope))
  face_to_bfs = Dict(F => bfs for (F, bfs) in b.bubbles)

  # Exponent β (|β| = r) of the lattice node x_β = Σᵢ βᵢ/r · Vᵢ carrying the
  # predof of a bubble entry.
  node_multiindex(::RotatingPΛBasis, k::Int, α) = copy(α)
  function node_multiindex(::TrimmedPΛBasis, e::Tuple{Int,Int}, α)
    β = copy(α); β[e[1]] += 1; return β
  end

  # Predof covector at x_β: the entry's own (physical, Cartesian) form
  # direction evaluated at the node.  λ_i(x_β) = β_i/r for the Whitney form.
  node_covector(b::RotatingPΛBasis, w::Int, β) = b.Ψ[w]
  function node_covector(b::TrimmedPΛBasis, w::Int, β)
    e1, e2 = b.e1s[w], b.e2s[w]
    return (β[e1]/r) * b.Je2[w] - (β[e2]/r) * b.Je1[w]
  end

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
    for (j, (w, k, α, _)) in enumerate(bfs)
      β = node_multiindex(b, k, α)
      rows[j] = get!(β_row, β) do
        push!(all_nodes, lattice_point(β))
        n = length(β_row) + 1; βs[n] = β; return n
      end
    end
    n_βs = length(β_row)
    moms = zeros(V, (n_βs, n_bubble))
    for (j, (w, k, α, _)) in enumerate(bfs)
      moms[rows[j], j] = node_covector(b, w, βs[rows[j]])
    end
    f_moments[gf]  = moms
    f_nodes[gf]    = (n_nodes+1) : (n_nodes+n_βs)
    f_own_moms[gf] = collect((n_dofs+1) : (n_dofs+n_bubble))
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

function _pλ_reffe(name::ReferenceFEName, ::Val{D}, r::Int) where D
  @assert D in (2, 3) "only D = 2, 3 supported, got D = $D"
  @assert r ≥ 1 "r must be ≥ 1, got $r"
  basis     = _pλ_basis(name, Val(D), r)
  polytope  = simplex_polytope(Val(D))
  n_dofs    = length(basis)
  predofs   = _pλ_predofs(basis, polytope, Val(D))
  face_dofs = get_face_own_funs(basis, polytope, CurlConformity())
  GenericRefFE{typeof(name)}(
    n_dofs, polytope, predofs, CurlConformity(), nothing, face_dofs, basis
  )
end

"""
    RotatingPΛRefFE(D, r) → GenericRefFE{RotatingPΛName,D}

H(curl)-conforming reference FE for the full rotating P_rΛ¹ basis on the
reference D-simplex (D = 2, 3; r ≥ 1).
"""
RotatingPΛRefFE(D::Int, r::Int) = _pλ_reffe(rotating_pλ, Val(D), r)

"""
    TrimmedPΛRefFE(D, r) → GenericRefFE{TrimmedPΛName,D}

H(curl)-conforming reference FE for the trimmed P_r⁻Λ¹ basis on the
reference D-simplex (D = 2, 3; r ≥ 1).
"""
TrimmedPΛRefFE(D::Int, r::Int) = _pλ_reffe(trimmed_pλ, Val(D), r)

# Standard factories: ReferenceFE(TRI, rotating_pλ, r), etc.
ReferenceFE(p::Polytope, ::RotatingPΛName, r::Int) = RotatingPΛRefFE(num_dims(p), r)
ReferenceFE(p::Polytope, ::TrimmedPΛName, r::Int) = TrimmedPΛRefFE(num_dims(p), r)

# ─────────────────────────────────────────────────────────────────────────────
# Geometric decomposition API (GeometricDecompositions.jl)
# ─────────────────────────────────────────────────────────────────────────────

has_geometric_decomposition(b::RotatingPΛBasis{D}, p::Polytope, conf::Conformity) where D =
  _pλ_geo_decomposition(b, Val(D), p, conf)
has_geometric_decomposition(b::TrimmedPΛBasis{D}, p::Polytope, conf::Conformity) where D =
  _pλ_geo_decomposition(b, Val(D), p, conf)

# Disambiguate against the generic (shapefuns, p, ::L2Conformity) method.
has_geometric_decomposition(::RotatingPΛBasis, ::Polytope, ::L2Conformity) = true
has_geometric_decomposition(::TrimmedPΛBasis,  ::Polytope, ::L2Conformity) = true

function _pλ_geo_decomposition(b, ::Val{D}, p, conf) where D
  conf isa L2Conformity && return true
  (!is_simplex(p) || D != num_dims(p)) && return false
  _are_barycoords_relative_to_simplex(b.scalar_bernstein_basis, p) ||
    return false
  conf isa CurlConformity
end

function get_face_own_funs(b::_PΛBases, p::Polytope, conf::CurlConformity)
  @assert has_geometric_decomposition(b, p, conf)
  faces = get_faces(p)
  face_own_funs = [Int[] for _ in 1:length(faces)]
  for (F, bubble_functions) in b.bubbles
    face = findfirst(face -> F ⊆ face, faces)
    w_first = first(bubble_functions)[1]
    w_last  = last(bubble_functions)[1]
    face_own_funs[face] = collect(w_first:w_last)
  end
  face_own_funs
end
