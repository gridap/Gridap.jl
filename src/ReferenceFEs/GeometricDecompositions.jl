################################
# Geometric decomposition APIs #
################################

"""
    has_geometric_decomposition(shapefuns, p::Polytope, ::Conformity) -> Bool

Tells whether `shapefuns` is a geometrically decomposed basis on `p` for the
given conformity. This is always true for `L2Conformity()`.

Otherwise, the decomposition is defined relatively to the appropriate trace:
- For `GradConformity()`, the trace is the restriction to the face, defined on all boundary faces,
- For `CurlConformity()`, the trace is the tangential trace to edges and tangential component to 2D facets,
- For `DivConformity()`, the trace is the normal trace to facets (rotated tangents in 2D).
"""
has_geometric_decomposition(shapefuns, p::Polytope, ::Conformity) = false
has_geometric_decomposition(shapefuns, p::Polytope, ::L2Conformity) = true

"""
    get_face_own_funs(shapefuns, p::Polytope, ::Conformity) -> Vector{Vector{Int}}

Essentially the same as [`get_face_own_dofs`](@ref), but for the `shapefuns`
basis instead of a DoF basis.

`shapefun` must implement the geometric decomposition on `p` for the given
conformity, this can be checked using [`has_geometric_decomposition`](@ref).
"""
get_face_own_funs(shapefuns, ::Polytope, ::Conformity) = @abstractmethod

function get_face_own_funs(shapefuns, p::Polytope, ::L2Conformity)
  _l2_conforming_own_funs(shapefuns,p)
end

function _l2_conforming_own_funs(shapefuns,p)
  r = [Int[] for _ in 1:num_faces(p)]
  r[end] = collect(1:length(shapefuns))
  r
end

"""
    get_facet_flux_sign_flip(shapefun, p::Polytope, ::DivConformity)

Return the (diagonal) change of basis matrix to make the flux of facet-owned
polynomials of `b` be oriented outwards the facet.

`shapefun` must implement the geometric decomposition on `p` for `DivConformity()`,
this can be checked using [`has_geometric_decomposition`](@ref).

# Extended help

The gluing of div-conforming shape functions assumes that all facet-owned
shapefuns have consistent orientation between facets (if the first shapefun of
facet 1 has outwards flux, all first shapefun of other facets must also have
outwards flux), see also [`NormalSignMap`](@ref Gridap.FESpaces.NormalSignMap).

This is not the case for the `BarycentricP(m)ΛBases` by default, their flux is
oriented like the sign of the permutation of the facet node indices.
"""
get_facet_flux_sign_flip(shapefuns, ::Polytope, ::L2Conformity) = @abstractmethod


##############################################################
# Geometric decompositions of barycentric bases on simplices #
##############################################################

# BernsteinBasisOnSimplex

function has_geometric_decomposition(
  b::BernsteinBasisOnSimplex{D}, p::Polytope, conf::Conformity) where D

  conf isa L2Conformity && return true

  !is_simplex(p) || D != num_dims(p) && return false
  if !_are_barycoords_relative_to_simplex(b, p)
    @warn """
      The barycentric coordinates of the given basis is not defined relative to the given simplex vertices.
      $(sprint(Base.show_backtrace, stacktrace()))
    """
    return false
  end

  conf isa H1Conformity && return true

  false
end

function get_face_own_funs(
  b::BernsteinBasisOnSimplex{D,V}, p::Polytope, conf::Conformity) where {D,V}

  @check has_geometric_decomposition(b,p,conf)

  conf isa L2Conformity && return _l2_conforming_own_funs(b,p)

  faces = get_faces(p)
  num_faces = length(faces)
  face_own_funs = Vector{Int}[ Int[] for _ in 1:num_faces]

  K = get_order(b)
  ncomp = num_indep_components(V)
  id = 1
  for α in bernstein_terms(K,D)
    F = findall(>(0), α)                    # vertices of the face owning x_α
    face = findfirst(face -> F⊆face, faces) # p face number
    # this should be guaranteed by has_geometric_decomposition
    isnothing(face) && @unreachable
    # the ncomp components holding the Bα scalar shapefun are contiguous due
    # to Polynomials._cartprod_set_value!
    append!(face_own_funs[face], id:id+ncomp-1)
    id += ncomp
  end

  face_own_funs
end

function _are_barycoords_relative_to_simplex(
  b::BernsteinBasisOnSimplex{D}, simplex::Polytope) where D

  @check is_simplex(simplex) && D == num_dims(simplex)

  vertices = get_vertex_coordinates(simplex)
  M = b.cart_to_bary_matrix
  vλ = [ Polynomials._cart_to_bary(v, M) for v in vertices ] # return SVectors ...

  T = eltype(eltype(vλ))
  V = VectorValue{D+1,T}
  I_cols = component_basis(V)
  vλ = reinterpret(V, vλ)

  vλ ≈ I_cols
end

# BarycentricP(m)ΛBasis

_BaryPΛBasis = Polynomials._BaryPΛBasis
function _is_rotated_90(ids::Polynomials.BarycentricPΛIndices)
  comps = ids.components
  length(comps) !== 2 && return false # only for 2D
  comps[2][3] < 0 # sign of second component flipped
end

function has_geometric_decomposition(b::_BaryPΛBasis, p::Polytope, conf::Conformity)
  D = get_dimension(b)
  k = b.k # form order

  conf isa L2Conformity && return true

  !is_simplex(p) || D != num_dims(p) && return false
  if !_are_barycoords_relative_to_simplex(b.scalar_bernstein_basis, p)
    @warn """
      The barycentric coordinates of the given basis is not defined relative to the given simplex vertices.
      $(sprint(Base.show_backtrace, stacktrace()))
    """
    return false
  end

  conf isa H1Conformity   && k == 0                     && return true
  is_rotated_90 = _is_rotated_90(b._indices)
  conf isa CurlConformity && k == 1   && !is_rotated_90 && return true
  correct_rotate = D==2 ? is_rotated_90 : true
  conf isa DivConformity  && k == D-1 && correct_rotate && return true

  false
end

function get_face_own_funs(b::_BaryPΛBasis, p::Polytope, conf::Conformity)
  @check has_geometric_decomposition(b,p,conf)

  conf isa L2Conformity && return _l2_conforming_own_funs(b,p)

  faces = get_faces(p)
  num_faces = length(faces)
  face_own_funs = Vector{Int}[ Int[] for _ in 1:num_faces]

  for (F, bubble_functions) in get_bubbles(b)
    face = findfirst(face -> F⊆face, faces)
    # this should be guaranteed by has_geometric_decomposition
    isnothing(face) && @unreachable
    # Polynomials._check_PΛ_indices guaranties bubble shapefun ids are contiguous
    w_first = first(bubble_functions)[1]
    w_last  = last( bubble_functions)[1]
    face_own_funs[face] = collect(w_first:w_last)
  end

  face_own_funs
end

function get_facet_flux_sign_flip(
  b::_BaryPΛBasis, p::Polytope{D}, conf::DivConformity) where D

  facet_range = get_dimrange(p,D-1)
  face_own_funs = get_face_own_funs(b,p,conf)
  sign_flip = MVector(tfill(1, Val(length(b)))...)

  for (face, own_funs) in enumerate(face_own_funs)
    if face ∈ facet_range
      sign_flip[own_funs] .= iseven(face-first(facet_range)) ? -1 : 1
      # Equivalent definition:
      # F = get_faces(p)[face][1:D]
      # sign_flip[own_funs] .= -Polynomials._combination_sign(F))
    end
  end

  sign_flip = Diagonal(sign_flip)
end
