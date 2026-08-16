# RotatingPLambda/PΛRefFEs.jl
#
# Reference FEs for the rotating P_rΛ¹ and trimmed P_r⁻Λ¹ bases: thin Gridap
# ReferenceFEs whose cross-cell H(curl) conformity is obtained through a
# per-cell change of basis given by the closed-form vertex-relabelling
# calculus of Gridap.Polynomials (RotatingPLambda/PΛRotations.jl and
# RotatingPLambda/PΛTrimmedRotations.jl). The per-cell change of basis itself
# (the compute_cell_bases_changes methods) lives in
# src/FESpaces/RotatingPLambda/PΛCellBasesChanges.jl.
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

using Gridap.Polynomials: RotatingPΛBasis, TrimmedPΛBasis, _PΛBases

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
# Reference-simplex helpers
# ─────────────────────────────────────────────────────────────────────────────

_pλ_polytope(::Val{2}) = TRI
_pλ_polytope(::Val{3}) = TET

# V_1 = origin, V_{j+1} = e_j (barycentric λ = (1−Σx, x…)).
_pλ_reference_vertices(::Val{D}) where D =
  ntuple(i -> Point(ntuple(j -> i == j+1 ? 1.0 : 0.0, D)), D+1)

# Sorted vertex ids of every Gridap face of the polytope, indexed by gf.
function _pλ_face_vertices(polytope, ::Val{D}) where D
  dimranges = get_dimranges(polytope)
  fv = Vector{Vector{Int}}(undef, num_faces(polytope))
  for d in 0:D
    for (li, verts) in enumerate(get_faces(polytope, d, 0))
      fv[dimranges[d+1][li]] = sort(collect(Int, verts))
    end
  end
  fv
end

_pλ_value_type(b::RotatingPΛBasis) = eltype(b.Ψ)
_pλ_value_type(b::TrimmedPΛBasis)  = eltype(b.Je1)

# Exponent β (|β| = r) of the BB lattice node x_β = Σᵢ βᵢ/r · Vᵢ carrying the
# predof of a bubble entry.
_pλ_node_multiindex(::RotatingPΛBasis, k::Int, α) = copy(α)
function _pλ_node_multiindex(::TrimmedPΛBasis, e::Tuple{Int,Int}, α)
  β = copy(α)
  β[e[1]] += 1
  β
end

# Predof covector at x_β: the entry's own (physical, Cartesian) form direction
# evaluated at the node. λ_i(x_β) = β_i/r for the trimmed Whitney form.
_pλ_node_covector(b::RotatingPΛBasis, w::Int, β, r) = b.Ψ[w]
function _pλ_node_covector(b::TrimmedPΛBasis, w::Int, β, r)
  e1, e2 = b.e1s[w], b.e2s[w]
  (β[e1]/r) * b.Je2[w] - (β[e2]/r) * b.Je1[w]
end

# ─────────────────────────────────────────────────────────────────────────────
# Predof data (MomentBasedDofBasis arrays), driven by the bubble tables
# ─────────────────────────────────────────────────────────────────────────────

function _pλ_predof_data(b::_PΛBases, polytope, ::Val{D}) where D
  r      = get_order(b)
  V      = _pλ_value_type(b)
  n_gf   = num_faces(polytope)
  fverts = _pλ_face_vertices(polytope, Val(D))
  verts  = _pλ_reference_vertices(Val(D))
  face_to_bfs = Dict(F => bfs for (F, bfs) in b.bubbles)

  all_nodes  = Point{D,Float64}[]
  f_moments  = [zeros(V, 0, 0) for _ in 1:n_gf]
  f_nodes    = [1:0 for _ in 1:n_gf]
  f_own_moms = [Int[] for _ in 1:n_gf]
  node_count = 0
  dof_count  = 0

  for gf in 1:n_gf
    bfs = get(face_to_bfs, fverts[gf], nothing)
    bfs === nothing && continue
    nw    = length(bfs)
    βs    = Vector{Vector{Int}}()
    β_row = Dict{Vector{Int},Int}()
    rows  = Vector{Int}(undef, nw)
    for (j, (w, k, α, _)) in enumerate(bfs)
      β = _pλ_node_multiindex(b, k, α)
      rows[j] = get!(β_row, β) do
        push!(βs, β)
        length(βs)
      end
    end
    moms = zeros(V, length(βs), nw)
    for (j, (w, k, α, _)) in enumerate(bfs)
      moms[rows[j], j] = _pλ_node_covector(b, w, βs[rows[j]], r)
    end
    for β in βs
      push!(all_nodes,
        Point(ntuple(a -> sum(β[i] * verts[i][a] for i in 1:D+1) / r, D)))
    end
    f_moments[gf]  = moms
    f_nodes[gf]    = node_count+1 : node_count+length(βs)
    f_own_moms[gf] = collect(dof_count+1 : dof_count+nw)
    node_count += length(βs)
    dof_count  += nw
  end
  @assert dof_count == length(b)

  all_nodes, f_moments, f_nodes, f_own_moms
end

# face_own_dofs, indexed by basis function id w (contiguous per face).
function _pλ_face_own_dofs(b::_PΛBases, polytope, ::Val{D}) where D
  fverts = _pλ_face_vertices(polytope, Val(D))
  face_to_bfs = Dict(F => bfs for (F, bfs) in b.bubbles)
  map(1:num_faces(polytope)) do gf
    bfs = get(face_to_bfs, fverts[gf], nothing)
    bfs === nothing && return Int[]
    ws = Int[bf[1] for bf in bfs]
    @assert ws == collect(first(ws):last(ws)) "bubble ids not contiguous on face"
    ws
  end
end

# Identity face-own-dof permutation for EVERY vertex permutation of every face:
# the rotation change of basis absorbs all orientation effects, so the
# pindex machinery must not reorder anything.
function _pλ_identity_dof_perms(polytope, face_own_dofs)
  vtx_perms = get_face_vertex_permutations(polytope)
  [[collect(1:length(own)) for _ in vtx_perms[gf]]
   for (gf, own) in enumerate(face_own_dofs)]
end

# ─────────────────────────────────────────────────────────────────────────────
# PΛRefFE: GenericRefFE wrapper carrying the identity permutation table
# ─────────────────────────────────────────────────────────────────────────────

"""
    PΛRefFE{Name,D} <: ReferenceFE{D}

Reference FE for the rotating P_rΛ¹ (`Name = RotatingPΛName`) or trimmed
P_r⁻Λ¹ (`Name = TrimmedPΛName`) basis on the reference D-simplex. H(curl)
conforming; shape functions ARE the basis (prebasis == shapefuns), the dof
basis is its dual (computed from pointwise BB-lattice moments). Wraps a
`GenericRefFE` to provide identity face-own-dof permutations for every
pindex — cross-cell conformity is delegated entirely to
`compute_cell_bases_changes` (see `Gridap.FESpaces`), i.e. to the rotation
calculus.
"""
struct PΛRefFE{Name,D} <: ReferenceFE{D}
  reffe     :: GenericRefFE{Name,D}
  dof_perms :: Vector{Vector{Vector{Int}}}
end

get_name(::Type{<:PΛRefFE{Name}}) where Name                  = Name()
Conformity(rf::PΛRefFE)                                       = Conformity(rf.reffe)
get_polytope(rf::PΛRefFE)                                     = get_polytope(rf.reffe)
get_prebasis(rf::PΛRefFE)                                     = get_prebasis(rf.reffe)
get_dof_basis(rf::PΛRefFE)                                    = get_dof_basis(rf.reffe)
num_dofs(rf::PΛRefFE)                                         = num_dofs(rf.reffe)
get_face_own_dofs(rf::PΛRefFE, c::Conformity)                 = get_face_own_dofs(rf.reffe, c)
get_face_dofs(rf::PΛRefFE)                                    = get_face_dofs(rf.reffe)
get_shapefuns(rf::PΛRefFE)                                    = get_shapefuns(rf.reffe)
get_metadata(rf::PΛRefFE)                                     = get_metadata(rf.reffe)
get_face_own_dofs_permutations(rf::PΛRefFE, ::Conformity)     = rf.dof_perms

function _pλ_reffe(name::ReferenceFEName, ::Val{D}, r::Int) where D
  @assert D in (2, 3) "only D = 2, 3 supported, got D = $D"
  @assert r ≥ 1 "r must be ≥ 1, got $r"
  basis     = _pλ_basis(name, Val(D), r)
  polytope  = _pλ_polytope(Val(D))
  n         = length(basis)
  nodes, f_moments, f_nodes, f_own_moms = _pλ_predof_data(basis, polytope, Val(D))
  predofs   = MomentBasedDofBasis(nodes, f_moments, f_nodes, f_own_moms)
  face_dofs = _pλ_face_own_dofs(basis, polytope, Val(D))
  Name      = typeof(name)
  # predofs-constructor: shapefuns = basis, dofs = compute_dofs(predofs, basis)
  # (the dual basis of the shape functions; exact duality is what makes
  # interpolation reproduce fields of the space).
  reffe = GenericRefFE{Name}(n, polytope, predofs, CurlConformity(),
                             nothing, face_dofs, basis)
  PΛRefFE{Name,D}(reffe, _pλ_identity_dof_perms(polytope, face_dofs))
end

"""
    RotatingPΛRefFE(D, r) → PΛRefFE{RotatingPΛName,D}

H(curl)-conforming reference FE for the full rotating P_rΛ¹ basis on the
reference D-simplex (D = 2, 3; r ≥ 1).
"""
RotatingPΛRefFE(D::Int, r::Int) = _pλ_reffe(rotating_pλ, Val(D), r)

"""
    TrimmedPΛRefFE(D, r) → PΛRefFE{TrimmedPΛName,D}

H(curl)-conforming reference FE for the trimmed P_r⁻Λ¹ basis on the
reference D-simplex (D = 2, 3; r ≥ 1).
"""
TrimmedPΛRefFE(D::Int, r::Int) = _pλ_reffe(trimmed_pλ, Val(D), r)

# Standard factories: ReferenceFE(TRI, rotating_pλ, r), etc.
ReferenceFE(p::Polytope, ::RotatingPΛName, r::Int) =
  RotatingPΛRefFE(num_dims(p), r)
ReferenceFE(p::Polytope, ::TrimmedPΛName, r::Int) =
  TrimmedPΛRefFE(num_dims(p), r)

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
