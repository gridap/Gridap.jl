
# MomentBasedDofBasis

struct Moment <: Dof end

"""
    struct MomentBasedDofBasis{P,V} <: AbstractVector{Moment}

Implementation of a basis of discretized moment DoFs, where `P` is the type of
the quadrature nodes, and `V` the value type of the shape functions.
"""
struct MomentBasedDofBasis{P,V} <: AbstractVector{Moment}
  nodes::Vector{P}
  face_moments::Vector{Array{V}}
  face_nodes::Vector{UnitRange{Int}}
  face_own_moms::Vector{Vector{Int}}

  function MomentBasedDofBasis(nodes,f_moments,f_nodes,f_own_moms)
    P = eltype(nodes)
    V = eltype(eltype(f_moments))
    new{P,V}(nodes,f_moments,f_nodes,f_own_moms)
  end

  # Unused and untested
  #function MomentBasedDofBasis(f_nodes,f_moments,f_own_moms)
  #  P = eltype(eltype(f_nodes))
  #  V = eltype(eltype(f_moments))
  #  nodes = P[]
  #  face_nodes = UnitRange{Int}[]
  #  nfaces = length(f_nodes)
  #  n = 1
  #  for fi in 1:nfaces
  #    nodes_fi = f_nodes[fi]
  #    nini = n
  #    for node_fi in nodes_fi
  #      push!(nodes,node_fi)
  #      n += 1
  #    end
  #    nend = n-1
  #    push!(face_nodes,nini:nend)
  #  end
  #  new{P,V}(nodes,f_moments,face_nodes,f_own_moms)
  #end
end

Base.size(a::MomentBasedDofBasis) = (num_dofs(a),)
Base.axes(a::MomentBasedDofBasis) = (Base.OneTo(num_dofs(a)),)
Base.getindex(a::MomentBasedDofBasis,i::Integer) = Moment()
Base.IndexStyle(::MomentBasedDofBasis) = IndexLinear()

"""
    get_nodes(b::MomentBasedDofBasis)

Get the vector of DoF quadrature nodes of `b`.
"""
get_nodes(b::MomentBasedDofBasis) = b.nodes

"""
    get_face_moments(b::MomentBasedDofBasis)

Return the vector of discretized moments for each face of the underlying polytope.
"""
get_face_moments(b::MomentBasedDofBasis) = b.face_moments

"""
    get_face_own_dofs(b::MomentBasedDofBasis)

The ownership of `b`'s dofs to faces of the underlying polytope is defined as
the face the moment was defined on.
"""
get_face_own_dofs(b::MomentBasedDofBasis) = b.face_own_moms

"""
    get_face_nodes_dofs(b::MomentBasedDofBasis)

Return the moment quadrature node indices on each face of the underlying polytope.
"""
get_face_nodes_dofs(b::MomentBasedDofBasis) = b.face_nodes

function num_dofs(b::MomentBasedDofBasis)
  n = 0
  for m in b.face_moments
    n += size(m,2)
  end
  n
end

function return_cache(b::MomentBasedDofBasis{P,V}, field) where {P,V}
  cf = return_cache(field,b.nodes)
  vals = evaluate!(cf,field,b.nodes)
  Vr = eltype(vals)
  T = typeof( zero(V) ⊙ zero(Vr) )
  r = Array{T}(undef, (num_dofs(b), size(field)...))
  c = CachedArray(r)
  return c, cf
end

function evaluate!(cache, b::MomentBasedDofBasis, field::Field)
  c, cf = cache
  setsize!(c, size(b))
  vals = evaluate!(cf,field,b.nodes)
  dofs = c.array

  o = 1
  z = zero(eltype(dofs))
  face_nodes = b.face_nodes
  face_moments = b.face_moments
  for face in eachindex(face_moments)
    moments = face_moments[face]
    if !iszero(length(moments))
      nodes = face_nodes[face]
      ni,nj = size(moments)
      for j in 1:nj
        dofs[o] = z
        for i in 1:ni
          dofs[o] += moments[i,j] ⊙ vals[nodes[i]]
        end
        o += 1
      end
    end
  end

  return dofs
end

function evaluate!(cache, b::MomentBasedDofBasis, field::AbstractVector{<:Field})
  c, cf = cache
  setsize!(c, (size(b,1),length(field)))
  vals = evaluate!(cf,field,b.nodes)
  dofs = c.array

  o = 1
  na = size(vals,2)
  z = zero(eltype(dofs))
  face_nodes = b.face_nodes
  face_moments = b.face_moments
  for face in eachindex(face_moments)
    moments = face_moments[face]
    if !iszero(length(moments))
      nodes = face_nodes[face]
      ni,nj = size(moments)
      for j in 1:nj
        for a in 1:na
          dofs[o,a] = z
          for i in 1:ni
            dofs[o,a] += moments[i,j] ⊙ vals[nodes[i],a]
          end
        end
        o += 1
      end
    end
  end

  return dofs
end

"""
    MomentBasedDofBasis(p::Polytope, prebasis::AbstractVector{<:Field}, moments)

Creates a basis of DoFs defined by moments on faces of `p`.

`moments` is a vector of moment descriptors, each one is given by a triplet
(f,σ,μ) where
  - f is collection of ids of faces Fₖ of `p`, that index `get_faces(p)`,
  - σ is a function σ(φ,μ,ds) **linear** in φ and μ that takes two Field-vectors φ and μ and a `FaceMeasure` ds and returns a Field-like object to be integrated over each face Fₖ,
  - μ is a polynomials basis on Fₖ.

The moment DoFs are thus defined by φ -> ∫_Fₖ σ(φ,μᵢ,ds)dFₖ,  ∀ σ,k,i.
In the final basis, DoFs are ordered by moment, then by face, then by "test" polynomial.

All the faces in a moment must be of the same type (have same reference face).
"""
function MomentBasedDofBasis(
    p::Polytope{D},
    prebasis::AbstractVector{<:Field},
    moments::AbstractVector{<:Tuple},
  ) where D

  n_faces = num_faces(p)
  n_moments = length(moments)
  face_dims = get_facedims(p)
  face_offsets = get_offsets(p)
  reffaces, face_types = _compute_reffaces_and_face_types(p)

  V = return_type(prebasis)
  φ_vec = representatives_of_componentbasis_dual(V)
  φ = map(constant_field,φ_vec)

  # Create face measures for each moment
  order = get_order(prebasis)
  measures = Vector{FaceMeasure}(undef,n_moments)
  for (k,(faces,σ,μ)) in enumerate(moments)
    ftype = face_types[first(faces)]
    @check all(isequal(ftype), face_types[faces])
    qdegree = order + get_order(μ) + 1
    fp = reffaces[ftype]
    measures[k] = FaceMeasure(p,fp,qdegree)
  end

  # Count number of moments and quad pts per face
  face_n_moms = zeros(Int,n_faces)
  face_n_nodes = zeros(Int,n_faces)
  for ((faces,σ,μ),ds) in zip(moments,measures)
    face_n_moms[faces] .+= length(μ)
    face_n_nodes[faces] .+= num_points(ds.quad)
  end

  # Compute face moment and node indices
  n_moms = 0
  n_nodes = 0
  face_nodes = Vector{UnitRange{Int}}(undef, n_faces)
  face_moments = Vector{Array{V}}(undef, n_faces)
  face_own_moms = Vector{Vector{Int}}(undef,n_faces)
  for face in 1:n_faces
    n_moms_i = face_n_moms[face]
    n_nodes_i = face_n_nodes[face]
    face_nodes[face] = (n_nodes+1):(n_nodes+n_nodes_i)
    face_moments[face] = zeros(V,n_nodes_i,n_moms_i)
    face_own_moms[face] = collect((n_moms+1):(n_moms+n_moms_i))
    n_moms += n_moms_i
    n_nodes += n_nodes_i
  end

  # Compute face moments and nodes
  fill!(face_n_moms,0)
  fill!(face_n_nodes,0)
  nodes = Vector{Point{D,Float64}}(undef,n_nodes)
  for ((faces,σ,μ),ds) in zip(moments,measures)
    cache = return_cache(σ,φ,μ,ds)

    for face in faces
      d = face_dims[face]
      lface = face - face_offsets[d+1]
      set_face!(ds,lface)

      # vals : (nN, nμ, nφ), coords : (nN)
      vals, coords = evaluate!(cache,σ,φ,μ,ds)
      # test_moment(σ,prebasis,μ,ds)

      mom_offset = face_n_moms[face]
      node_offset = first(face_nodes[face]) + face_n_nodes[face] - 1
      for i in axes(vals,1)
        for j in axes(vals,2)
          face_moments[face][i,j+mom_offset] = V(vals[i,j,:]...)
        end
        nodes[i+node_offset] = coords[i]
      end

      face_n_nodes[face] += size(vals,1)
      face_n_moms[face] += size(vals,2)
    end
  end

  MomentBasedDofBasis(nodes, face_moments, face_nodes, face_own_moms)
end

# Unused and untested
#function test_moment(σ,prebasis,μ,ds)
#  T = return_type(prebasis)
#  φ = map(constant_field,dual_component_basis_representatives(T))
#  vals, coords = evaluate(σ,φ,μ,ds)
#
#  φx = evaluate(prebasis, coords) # (nN, nφ)
#  σx, _ = evaluate(σ,prebasis,μ,ds)
#
#  σx_bis = zeros(size(vals,1),size(vals,2),size(σx,3))
#  for i in axes(vals,1)
#    for j in axes(vals,2)
#      cx = T(vals[i,j,:]...)
#      σx_bis[i,j,:] .= map(y -> inner(y,cx),φx[i,:])
#    end
#  end
#
#  println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
#
#  println(" > Values: ")
#  for i in axes(vals,2)
#    println("    >> Mu = ",i)
#    for j in axes(vals,3)
#      println("     >>> Phi -> ", sum(vals[:,i,j]))
#    end
#  end
#
#  println(" > Moments: ")
#  for i in axes(vals,1)
#    display(σx[i,:,:])
#    display(σx_bis[i,:,:])
#    println("----------------------------------------")
#  end
#  @assert σx ≈ σx_bis
#  println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
#end


###############
# FaceMeasure #
###############

# evaluate! and return_cache arn't type stable, at least because `quad` and
# `fmaps` are not concretely typed.
# I think we could use a MVector{1,Int} to store `face` and make the
# struct non-mutable, this would likely be more performent
mutable struct FaceMeasure{Df,Dc}
  face ::Int
  cpoly::Polytope{Dc}
  fpoly::Polytope{Df}
  quad ::Quadrature
  fmaps::Vector{<:Field}
  function FaceMeasure(
    cpoly::Polytope{Dc},fpoly::Polytope{Df},order::Int
  ) where {Df,Dc}
    # Quadrature on the face
    quad = Quadrature(fpoly,order)
    # Face to cell coordinate map
    #if Df == Dc
    #  fmaps = [GenericField(identity)]
    #else # TODO: Could this be an AffineMap?
      fcoords = get_face_coordinates(cpoly,Df)
      basis = get_shapefuns(LagrangianRefFE(Float64,fpoly,1))
      fmaps = map(c -> linear_combination(c,basis),fcoords)
    #end
    new{Df,Dc}(1,cpoly,fpoly,quad,fmaps)
  end
end

function set_face!(m::FaceMeasure{Df},face::Int) where {Df}
  @assert 0 < face <= num_faces(m.cpoly, Df)
  m.face = face
  return m
end

# TODO: Normals are accesed, but tangent are computed on demand. This means
# that we will be repeating work unless we cache them.
function get_facet_normal(m::FaceMeasure{Df,Dc}) where {Df,Dc}
  @assert Df == Dc - 1
  n = get_facet_normal(m.cpoly)
  return ConstantField(n[m.face])
end

function get_edge_tangent(m::FaceMeasure{1,Dc}) where {Dc}
  t = get_edge_tangent(m.cpoly)
  return ConstantField(t[m.face])
end

# Matrix of the contravariant piola map from `m.fpoly` to to the face `m.face`
# of `m.cpoly`, used to extend a Df-dimensional vector in `m.fpoly` to a
# Dc-dimensional one that lives in the tangent space of the Dc-embedded Df-dimensional
# manifold `m.face`.
function get_extension(m::FaceMeasure{Df,Dc}) where {Df,Dc}
  @assert Df == Dc - 1
  vs = ReferenceFEs._nfaces_vertices(Float64,m.cpoly,Df)[m.face]
  J = TensorValue(hcat([vs[2]-vs[1]...],[vs[3]-vs[1]...]))
  return ConstantField(J/meas(transpose(J)))
end
# function get_extension(m::FaceMeasure{Df,Dc}) where {Df,Dc}
#   @assert Df == Dc - 1
#   fmap = m.fmaps[m.face]
#   J = Broadcasting(∇)(fmap)
#   return Operation(*)(Operation(transpose)(J),Operation(x -> 1/meas(x))(J))
# end

function Arrays.return_cache(
  σ::Function,                # σ(φ,μ,ds) -> Field/Array{Field}
  φ::AbstractArray{<:Field},  # φ: prebasis (defined on the cell)
  μ::AbstractArray{<:Field},  # μ: polynomial basis (defined on the face)
  ds::FaceMeasure             # ds: face measure
)
  fmap = ds.fmaps[ds.face]
  φf = transpose(Broadcasting(Operation(∘))(φ,fmap))
  return_cache(σ(φf,μ,ds),ds)
end

function Arrays.evaluate!(
  cache,
  σ::Function,                # σ(φ,μ,ds) -> Field/Array{Field}
  φ::AbstractArray{<:Field},  # φ: prebasis (defined on the cell)
  μ::AbstractArray{<:Field},  # μ: polynomial basis (defined on the face)
  ds::FaceMeasure             # ds: face measure
)
  fmap = ds.fmaps[ds.face]
  φf = transpose(Broadcasting(Operation(∘))(φ,fmap))
  evaluate!(cache,σ(φf,μ,ds),ds)
end

function Arrays.return_cache(f,ds::FaceMeasure)
  fmap = ds.fmaps[ds.face]

  xf = get_coordinates(ds.quad)
  w = get_weights(ds.quad)
  fmap_cache = return_cache(fmap,xf)

  detJ = Broadcasting(Operation(meas))(Broadcasting(∇)(fmap))
  detJ_cache = return_cache(detJ,xf)

  f_cache = return_cache(f,xf)
  return fmap_cache, detJ_cache, f_cache, xf, w
end

function Arrays.evaluate!(cache,f,ds::FaceMeasure)
  fmap_cache, detJ_cache, f_cache, xf, w = cache
  fmap = ds.fmaps[ds.face]

  detJ = Broadcasting(Operation(meas))(Broadcasting(∇)(fmap))
  dF = evaluate!(detJ_cache,detJ,xf)

  xc = evaluate!(fmap_cache,fmap,xf) # quad pts on the cell
  fx = evaluate!(f_cache,f,xf) # f evaluated on the quad pts
  fx .= (w .* dF) .* fx
  return fx, xc
end


##########################
# MomentBasedReferenceFE #
##########################

"""
    MomentBasedReferenceFE(
      name::ReferenceFEName,
      p::Polytope,
      prebasis::AbstractVector{<:Field},
      moments::AbstractVector{<:Tuple},
      conformity::Conformity;
      sh_is_pb=false
    )

Constructs a ReferenceFEs on `p` with a moment DoF (pre-)basis defined by
`moments`. See [`MomentBasedDofBasis`](@ref) for the specification of `moments`.

If `sh_is_pb=true`, `prebasis` is used as shape functions.
This requires it to fullfill a geometric decomposition relative to the faces of
`p` for `conformity`, see the [*Geometric decompositions*](@ref "Geometric decompositions")
section in the docs.

Warning, this function does not check that the moments are properly defined to implement
the given `conformity`. It is assumed that if ``σ(φ) = 0`` for a moment ``σ``
owned by a face ``f`` of `p`, then the `conformity`-trace of ``φ`` over ``f`` is zero.
"""
function MomentBasedReferenceFE(
  name::ReferenceFEName,
  p::Polytope,
  prebasis::AbstractVector{<:Field},
  moments::AbstractVector{<:Tuple},
  conformity::Conformity;
  sh_is_pb=false
)

  dof_basis = MomentBasedDofBasis(p, prebasis, moments)
  n_dofs = length(dof_basis)
  metadata = nothing

  if sh_is_pb
    # Use geometric decomposition to avoid shapefuns change of basis.
    # The dofs are computed as dual to the shapefuns, using a change of basis
    # from the moment basis

    @assert has_geometric_decomposition(prebasis, p, conformity)
    shapefuns, predofs = prebasis, dof_basis

    if conformity isa DivConformity
      sign_flip = get_facet_flux_sign_flip(shapefuns, p, conformity)
      shapefuns = linear_combination(sign_flip,shapefuns)
    end

    dofs = compute_dofs(predofs, shapefuns) # TODO smart inverse
    face_own_dofs = get_face_own_funs(prebasis, p, conformity)
    return GenericRefFE{typeof(name)}(
      n_dofs, p, predofs, conformity, metadata, face_own_dofs, shapefuns, dofs
    )
  end

  # else, standard prebasis inversion
  face_own_dofs = get_face_own_dofs(dof_basis)
  GenericRefFE{typeof(name)}(
    n_dofs, p, prebasis, dof_basis, conformity, metadata, face_own_dofs
  )
end

# Default polynomial type for moment based reference FEs
function _mom_reffe_default_PT(p)
  is_simplex(p) && return Bernstein
  is_n_cube(p) && return Legendre
  Monomial
end

