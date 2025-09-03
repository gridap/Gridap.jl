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

  conf isa GradConformity && return true

  false
end

function get_face_own_funs(
  b::BernsteinBasisOnSimplex{D,V}, p::Polytope, conf::GradConformity) where {D,V}

  @check has_geometric_decomposition(b,p,conf)

  faces = get_faces(p)
  n_faces = length(faces)
  face_own_funs = Vector{Int}[ Int[] for _ in 1:n_faces]

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

  conf isa GradConformity && k == 0                     && return true
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
  n_faces = length(faces)
  face_own_funs = Vector{Int}[ Int[] for _ in 1:n_faces]

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


###############################################################
# Geometric decompositions of tensor product bases on n-cubes #
###############################################################

#Polynomial bases admitting a 1D geometric decomposition on the SEGMENT, currently `ModalC0` and `Bernstein`.
const GD_1D_PT = Union{Polynomials.ModalC0, Bernstein}

# says which poly of the Kth order 1D basis does the SEGMENT vertices own
_SEGMENT_vertex_own_fun(::Type{Polynomials.ModalC0}, K) = (1, 2) # those are 1-x and x
_SEGMENT_vertex_own_fun(::Type{Bernstein}, K) = (1, K+1)         # those are (1-x)ᴷ and xᴷ

const _V0 = 1 # SEGMENT first vertex
const _V1 = 2 # SEGMENT second vertex
const _Vi = 3 # SEGMENT interior

function _compute_fixed_coords_to_face(p,::Val{D}) where D
  face_vertices = get_face_coordinates(p)
  fixed_coords_to_face = Dict{NTuple{D,Int},Int}()
  for (face,face_verts) in enumerate(face_vertices)
    fixed_coords = ntuple( i ->
        all(v->iszero(v[i]),face_verts) ? _V0 :
        all(v->isone( v[i]),face_verts) ? _V1 :
                                          _Vi, Val(D))
    fixed_coords_to_face[fixed_coords] = face
  end
  fixed_coords_to_face
end

@inline function _is_poly_reference_D_cube(p,D)
  !(is_n_cube(p) && D == num_dims(p)) && return false
  if D<4
    DCUBE = (VERTEX, SEGMENT, QUAD, HEX)[D+1]
  else
    DCUBE = ExtrusionPolytope(tfill(HEX_AXIS,Val(D)))
  end
  # Our 1D polynomial evaluations are decomposed in [0,1]ᴰ, so the vertices of
  # p must be the same as the Reference D-cube
  DCUBE_vertices = get_vertex_coordinates(DCUBE)
  p_vertices = get_vertex_coordinates(p)
  all(∈(DCUBE_vertices), p_vertices)
end


# CartProdPolyBasis

function has_geometric_decomposition(
  b::CartProdPolyBasis{D,V,<:GD_1D_PT}, p::Polytope, conf::Conformity) where {D,V}

  conf isa L2Conformity && return true

  !_is_poly_reference_D_cube(p,D) && return false

  conf isa GradConformity && minimum(b.orders) > 0 && return true

  # # could be generalized to Curl and Div conformity in this case, although
  # # it is quite redundant with `CompWiseTensorPolyBasis` if terms filtering
  # # is implemented for it
  # V <: VectorValue{D}     && minimum(b.orders) > 0 && return true

  false
end

function get_face_own_funs(
  b::CartProdPolyBasis{D,V,PT}, p::Polytope, conf::GradConformity) where {D,V,PT<:GD_1D_PT}

  @check has_geometric_decomposition(b,p,conf)

  face_own_funs = Vector{Int}[ Int[] for _ in 1:num_faces(p) ]
  fixed_coords_to_face = _compute_fixed_coords_to_face(p,Val(D))

  K = get_order(b)
  s0_owned, s1_owned = _SEGMENT_vertex_own_fun(PT,K)
  ncomp = num_indep_components(V)
  id = 1
  for ci in b.terms
    own_coords = ntuple( i -> ci[i] == s0_owned ? _V0 : ci[i] == s1_owned ? _V1 : _Vi, Val(D))
    face = fixed_coords_to_face[own_coords]
    append!(face_own_funs[face], id:id+ncomp-1)
    id += ncomp
  end

  face_own_funs
end

# CompWiseTensorPolyBasis

function has_geometric_decomposition(
  b::CompWiseTensorPolyBasis{D,V,<:GD_1D_PT}, p::Polytope, conf::Conformity) where {D,V}

  conf isa L2Conformity && return true

  !(V <: VectorValue{D}) && return false
  !_is_poly_reference_D_cube(p,D) && return false

  orders = MMatrix{D,D}(b.orders)
  conf isa GradConformity && minimum(orders)   > 0 && return true
  conf isa CurlConformity && minimum(orders+I) > 0 && return true
  conf isa DivConformity  && minimum(diag(orders)) > 0 && return true

  false
end

function get_face_own_funs(
  b::CompWiseTensorPolyBasis{D,V,PT}, p::Polytope, conf::Conformity) where {D,V,PT<:GD_1D_PT}

  @check has_geometric_decomposition(b,p,conf)
  conf isa L2Conformity && return _l2_conforming_own_funs(b,p)

  face_own_funs = Vector{Int}[ Int[] for _ in 1:num_faces(p) ]
  fixed_coords_to_face = _compute_fixed_coords_to_face(p,Val(D))

  # For Curl and Div conformity, some faces cannot own any shape functions,
  # e.g. the faces orthogonal to eₓ cannot own a Curl-conform shape function with
  # non-zero eₓ components, so _compute_mask(Val(D), 1, CurlConformity()) returns
  #   (true, false, ..., false)
  # to indicate that the ownership along x-axis is ignored for the first component.
  # This also ensures that no vertex can own a Curl-conforming function and that
  # only facets and interior can own a Div-conforming function, as expected.
  function _compute_mask(VD, d, conf)
    conf isa GradConformity && return MVector(tfill(true, VD))
    conf isa CurlConformity && return MVector(ntuple(i -> i==d, VD))
    conf isa DivConformity  && return MVector(ntuple(i -> i!=d, VD))
  end

  K = get_order(b)
  s0_owned, s1_owned = _SEGMENT_vertex_own_fun(PT,K)
  id = 1
  for (d,terms) in enumerate(Polynomials.get_comp_terms(b))
    mask = _compute_mask(Val(D),d,conf)
    for ci in terms
      own_coords = ntuple(
        i ->           mask[i] ? _Vi :
             ci[i] == s0_owned ? _V0 :
             ci[i] == s1_owned ? _V1 :
                                 _Vi, Val(D))
      face = fixed_coords_to_face[own_coords]

      push!(face_own_funs[face], id)
      id += 1
    end
  end
  face_own_funs
end

function get_facet_flux_sign_flip(
  b::CompWiseTensorPolyBasis{D,V,PT}, p::Polytope, conf::Conformity) where {D,V,PT<:GD_1D_PT}

  facet_range = get_dimrange(p,D-1)
  face_own_funs = get_face_own_funs(b,p,conf)
  sign_flip = MVector(tfill(1, Val(length(b)))...)

  for (face, own_funs) in enumerate(face_own_funs)
    if face ∈ facet_range
      # empirically determined, see tests
      sign_flip[own_funs] .= iseven(face-first(facet_range)) ? -1 : 1
    end
  end

  sign_flip = Diagonal(sign_flip)
end


