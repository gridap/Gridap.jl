struct DivConformity <: Conformity end

"""
    RaviartThomasRefFE(::Type{et},p::Polytope,order::Integer) where et

The `order` argument has the following meaning: the divergence of the  functions in this basis
is in the Q space of degree `order`.
"""
function RaviartThomasRefFE(::Type{et},p::Polytope,order::Integer) where et

  D = num_dims(p)

  prebasis = QCurlGradMonomialBasis{D}(et,order)

  nf_nodes, nf_moments = _RT_nodes_and_moments(et,p,order)
  #nf_nodes, nf_moments = _RT_nodes_and_moments(et,p,order,GenericField(identity))

  face_own_dofs = _face_own_dofs_from_moments(nf_moments)

  face_dofs = face_own_dofs

  dof_basis = MomentBasedDofBasis(nf_nodes, nf_moments)

  ndofs = num_dofs(dof_basis)

  #metadata = (et,order)
  metadata = nothing

  reffe = GenericRefFE{:RaviartThomas}(
    ndofs,
    p,
    prebasis,
    dof_basis,
    DivConformity(),
    metadata,
    face_dofs)

  reffe
end

function ReferenceFE(p::Polytope,::Val{:RaviartThomas}, order)
  RaviartThomasRefFE(Float64,p,order)
end

function ReferenceFE(p::Polytope,::Val{:RaviartThomas},::Type{T}, order) where T
  RaviartThomasRefFE(T,p,order)
end

function Conformity(reffe::GenericRefFE{:RaviartThomas},sym::Symbol)
  hdiv = (:Hdiv,:HDiv)
  if sym == :L2
    L2Conformity()
  elseif sym in hdiv
    DivConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a Raviart Thomas reference FE.

    Possible values of conformity for this reference fe are $((:L2, hdiv...)).
    """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{:RaviartThomas}, conf::DivConformity)
  get_face_dofs(reffe)
end

# Make it work for quads for the moment
function get_dof_basis(reffe::GenericRefFE{:RaviartThomas},phi::Field)
  get_dof_basis(reffe)
end

## First implement this:
#function get_dof_basis(reffe::GenericRefFE{:RaviartThomas},phi::Field)
#   p = get_polytope(reffe)
#   et, order = reffe.metadata
#   nf_nodes, nf_moments = _RT_nodes_and_moments(et,p,order,phi)
#   MomentBasedDofBasis(nf_nodes, nf_moments)
#end

# Then optimize this:
#function get_dof_basis(reffe::GenericRefFE{:RaviartThomas},phi::Field)
#  cache = return_cache(get_dof_basis,reffe,phi)
#  evaluate!(cache,get_dof_basis,reffe,phi)
#end
#
#function return_cache(::typeof(get_dof_basis),reffe::GenericRefFE{:RaviartThomas},phi::Field)
#  @notimplemented
#end
#
#function evaluate!(cache,::typeof(get_dof_basis),reffe::GenericRefFE{:RaviartThomas},phi::Field)
#  @notimplemented
#end

function _RT_nodes_and_moments(::Type{et}, p::Polytope, order::Integer) where et

  @notimplementedif ! is_n_cube(p)

  D = num_dims(p)
  ft = VectorValue{D,et}
  pt = Point{D,et}

  nf_nodes = [ zeros(pt,0) for face in 1:num_faces(p)]
  nf_moments = [ zeros(ft,0,0) for face in 1:num_faces(p)]

  fcips, fmoments = _RT_face_values(p,et,order)
  frange = get_dimrange(p,D-1)
  nf_nodes[frange] = fcips
  nf_moments[frange] = fmoments

  if (order > 0)
    ccips, cmoments = _RT_cell_values(p,et,order)
    crange = get_dimrange(p,D)
    nf_nodes[crange] = ccips
    nf_moments[crange] = cmoments
  end

  nf_nodes, nf_moments
end

# Ref FE to faces geomaps
function _ref_face_to_faces_geomap(p,fp)
  cfvs = get_face_coordinates(p,num_dims(fp))
  nc = length(cfvs)
  freffe = LagrangianRefFE(Float64,fp,1)
  fshfs = get_shapefuns(freffe)
  cfshfs = fill(fshfs, nc)
  # fgeomap = lincomb(cfshfs,cfvs)
  fgeomap = lazy_map(linear_combination,cfvs,cfshfs)
end

function _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)
  nc = length(fgeomap)
  c_fips = fill(fips,nc)
  c_wips = fill(wips,nc)
  pquad = lazy_map(evaluate,fgeomap,c_fips)
  c_fips, pquad, c_wips
end

function _broadcast(::Type{T},n,b) where T
  c = Array{T}(undef,size(b))
  for (ii, i) in enumerate(b)
    c[ii] = i*n
  end
  return c
end

function _RT_face_moments(p, fshfs, c_fips, fcips, fwips)
  nc = length(c_fips)
  cfshfs = fill(fshfs, nc)
  cvals = lazy_map(evaluate,cfshfs,c_fips)
  cvals = [fwips[i].*cvals[i] for i in 1:nc]
  # fns, os = get_facet_normal(p)
  fns = get_facet_normal(p)
  os = get_facet_orientations(p)
  # @santiagobadia : Temporary hack for making it work for structured hex meshes
  cvals = [ _broadcast(typeof(n),n*o,b) for (n,o,b) in zip(fns,os,cvals)]
  #cvals = [ _broadcast(typeof(n),n,b) for (n,b) in zip(fns,cvals)]
  return cvals
end

# It provides for every face the nodes and the moments arrays
function _RT_face_values(p,et,order)

  # Reference facet
  @assert is_simplex(p) || is_n_cube(p) "We are assuming that all n-faces of the same n-dim are the same."
  fp = Polytope{num_dims(p)-1}(p,1)

  # geomap from ref face to polytope faces
  fgeomap = _ref_face_to_faces_geomap(p,fp)

  # Nodes are integration points (for exact integration)
  # Thus, we define the integration points in the reference
  # face polytope (fips and wips). Next, we consider the
  # n-face-wise arrays of nodes in fp (constant cell array c_fips)
  # the one of the points in the polytope after applying the geopmap
  # (fcips), and the weights for these nodes (fwips, a constant cell array)
  # Nodes (fcips)
  degree = (order+1)*2
  fquad = Quadrature(fp,degree)
  fips = get_coordinates(fquad)
  wips = get_weights(fquad)

  c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

  # Moments (fmoments)
  # The RT prebasis is expressed in terms of shape function
  fshfs = MonomialBasis(et,fp,order)

  # Face moments, i.e., M(Fi)_{ab} = q_RF^a(xgp_RFi^b) w_Fi^b n_Fi ⋅ ()
  fmoments = _RT_face_moments(p, fshfs, c_fips, fcips, fwips)

  return fcips, fmoments

end

function _RT_cell_moments(p, cbasis, ccips, cwips)
  # Interior DOFs-related basis evaluated at interior integration points
  ishfs_iips = evaluate(cbasis,ccips)
  return cwips.*ishfs_iips
end

# It provides for every cell the nodes and the moments arrays
function _RT_cell_values(p,et,order)
  # Compute integration points at interior
  degree = 2*(order+1)
  iquad = Quadrature(p,degree)
  ccips = get_coordinates(iquad)
  cwips = get_weights(iquad)

  # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
  cbasis = QGradMonomialBasis{num_dims(p)}(et,order-1)
  cmoments = _RT_cell_moments(p, cbasis, ccips, cwips )

  return [ccips], [cmoments]

end

function _face_own_dofs_from_moments(f_moments)
  face_dofs = Vector{Int}[]
  o = 1
  for moments in f_moments
    ndofs = size(moments,2)
    dofs = collect(o:(o+ndofs-1))
    push!(face_dofs,dofs)
    o += ndofs
  end
  face_dofs
end

struct Moment <: Dof end

struct MomentBasedDofBasis{P,V} <: AbstractVector{Moment}
  nodes::Vector{P}
  face_moments::Vector{Array{V}}
  face_nodes::Vector{UnitRange{Int}}

  function MomentBasedDofBasis(nodes,f_moments,f_nodes)
    P = eltype(nodes)
    V = eltype(eltype(f_moments))
    new{P,V}(nodes,f_moments,f_nodes)
  end

  function MomentBasedDofBasis(f_nodes,f_moments)
    P = eltype(eltype(f_nodes))
    V = eltype(eltype(f_moments))
    nodes = P[]
    face_nodes = UnitRange{Int}[]
    nfaces = length(f_nodes)
    n = 1
    for fi in 1:nfaces
      nodes_fi = f_nodes[fi]
      nini = n
      for node_fi in nodes_fi
        push!(nodes,node_fi)
        n += 1
      end
      nend = n-1
      push!(face_nodes,nini:nend)
    end
    new{P,V}(nodes,f_moments,face_nodes)
  end
end

@inline Base.size(a::MomentBasedDofBasis) = (length(a.nodes),)
@inline Base.axes(a::MomentBasedDofBasis) = (axes(a.nodes,1),)
# @santiagobadia : Not sure we want to create the moment dofs
@inline Base.getindex(a::MomentBasedDofBasis,i::Integer) = Moment()
@inline Base.IndexStyle(::MomentBasedDofBasis) = IndexLinear()

get_nodes(b::MomentBasedDofBasis) = b.nodes
get_face_moments(b::MomentBasedDofBasis) = b.face_moments
get_face_nodes_dofs(b::MomentBasedDofBasis) = b.face_nodes

function num_dofs(b::MomentBasedDofBasis)
  n = 0
  for m in b.face_moments
    n += size(m,2)
  end
  n
end

function return_cache(b::MomentBasedDofBasis,field)
  cf = return_cache(field,b.nodes)
  vals = evaluate!(cf,field,b.nodes)
  ndofs = num_dofs(b)
  r = _moment_dof_basis_cache(vals,ndofs)
  c = CachedArray(r)
  (c, cf)
end

function _moment_dof_basis_cache(vals::AbstractVector,ndofs)
  T = eltype(vals)
  r = zeros(eltype(T),ndofs)
end

function _moment_dof_basis_cache(vals::AbstractMatrix,ndofs)
  _, npdofs = size(vals)
  T = eltype(vals)
  r = zeros(eltype(T),ndofs,npdofs)
end

function evaluate!(cache,b::MomentBasedDofBasis,field)
  c, cf = cache
  vals = evaluate!(cf,field,b.nodes)
  dofs = c.array
  _eval_moment_dof_basis!(dofs,vals,b)
  dofs
end

function _eval_moment_dof_basis!(dofs,vals::AbstractVector,b)
  o = 1
  z = zero(eltype(dofs))
  face_nodes = b.face_nodes
  face_moments = b.face_moments
  for face in 1:length(face_moments)
    moments = face_moments[face]
    if length(moments) != 0
      nodes = face_nodes[face]
      ni,nj = size(moments)
      for j in 1:nj
        dofs[o] = z
        for i in 1:ni
          dofs[o] += moments[i,j]⋅vals[nodes[i]]
        end
        o += 1
      end
    end
  end
end

function _eval_moment_dof_basis!(dofs,vals::AbstractMatrix,b)
  o = 1
  na = size(vals,2)
  z = zero(eltype(dofs))
  face_nodes = b.face_nodes
  face_moments = b.face_moments
  for face in 1:length(face_moments)
    moments = face_moments[face]
    if length(moments) != 0
      nodes = face_nodes[face]
      ni,nj = size(moments)
      for j in 1:nj
        for a in 1:na
          dofs[o,a] = z
          for i in 1:ni
            dofs[o,a] += moments[i,j]⋅vals[nodes[i],a]
          end
        end
        o += 1
      end
    end
  end
end

struct ContraVariantPiolaMap <: PushForwardMap end

function evaluate!(cache,::ContraVariantPiolaMap,v::Number,J::Number,detJ::Number)
  (1/detJ)*J⋅v
end

function evaluate!(cache,k::ContraVariantPiolaMap,v::AbstractVector{<:Field},phi::Field)
  Jt = ∇(phi)
  detJ = Operation(det)(Jt)
  J = Operation(transpose)(Jt)
  Broadcasting(Operation(k))(v,J,detJ)
end

#function lazy_map(
#  ::typeof(evaluate),
#  cell_shapefuns::LazyArray{<:Fill{<:ContraVariantPiolaMap}},
#  cell_qs::AbstractArray)
#
#  cell_ref_shapefuns = cell_shapefuns.f[1]
#  cell_map = cell_shapefuns.f[2]
#
#  cell_ref_shapefuns_q = lazy_map(evaluate,cell_ref_shapefuns,cell_qs)
#  cell_Jt = lazy_map(∇,cell_map)
#  cell_Jt_q = lazy_map(evaluate,cell_Jt,cell_qs)
#  cell_detJ_q = lazy_map(Broadcasting(det),cell_Jt_q)
#  cell_J_q = lazy_map(Broadcasting(transpose),cell_Jt_q)
#
#  k = cell_shapefuns.g.value
#  lazy_map(Broadcasting(k),cell_ref_shapefuns_q,cell_J_q,cell_detJ_q)
#end

function lazy_map(
  k::ContraVariantPiolaMap,
  cell_ref_shapefuns::AbstractArray{<:AbstractArray{<:Field}},
  cell_map::AbstractArray{<:Field})

  cell_Jt = lazy_map(∇,cell_map)
  cell_detJ = lazy_map(Operation(det),cell_Jt)
  cell_J = lazy_map(Operation(transpose),cell_Jt)

  lazy_map(Broadcasting(Operation(k)),cell_ref_shapefuns,cell_J,cell_detJ)
end

function evaluate!(cache,::ContraVariantPiolaMap,s::MomentBasedDofBasis,phi::Field)
  @notimplemented
  # phi_q = evaluate(phi,s.nodes)
  #moments = #from s.moments and phi_q
  #MomentBasedDofBasis(s.nodes,moments)
  ## More obvious but not so efficient:
  ##nf_nodes, nf_moments = _RT_nodes_and_moments(et,p,order,phi)
  ##dof_basis = MomentBasedDofBasis(nf_nodes, nf_moments,GenericField(identity))
end

PushForwardMap(reffe::GenericRefFE{:RaviartThomas}) = ContraVariantPiolaMap()



