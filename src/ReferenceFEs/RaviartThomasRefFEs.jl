struct DivConformity <: Conformity end

struct RaviartThomas <: ReferenceFEName end

const raviart_thomas = RaviartThomas()

"""
    RaviartThomasRefFE(::Type{et},p::Polytope,order::Integer) where et

The `order` argument has the following meaning: the divergence of the  functions in this basis
is in the Q space of degree `order`.
"""
function RaviartThomasRefFE(::Type{et},p::Polytope,order::Integer) where et

  D = num_dims(p)

  if is_n_cube(p)
    prebasis = QCurlGradMonomialBasis{D}(et,order)
  elseif is_simplex(p)
    prebasis = PCurlGradMonomialBasis{D}(et,order)
  else
    @notimplemented "H(div) Reference FE only available for cubes and simplices"
  end

  nf_nodes, nf_moments = _RT_nodes_and_moments(et,p,order,GenericField(identity))

  face_own_dofs = _face_own_dofs_from_moments(nf_moments)

  face_dofs = face_own_dofs

  dof_basis = MomentBasedDofBasis(nf_nodes, nf_moments)

  ndofs = num_dofs(dof_basis)

  metadata = nothing

  reffe = GenericRefFE{RaviartThomas}(
    ndofs,
    p,
    prebasis,
    dof_basis,
    DivConformity(),
    metadata,
    face_dofs)

  reffe
end

function ReferenceFE(p::Polytope,::RaviartThomas, order)
  RaviartThomasRefFE(Float64,p,order)
end

function ReferenceFE(p::Polytope,::RaviartThomas,::Type{T}, order) where T
  RaviartThomasRefFE(T,p,order)
end

function Conformity(reffe::GenericRefFE{RaviartThomas},sym::Symbol)
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

function get_face_own_dofs(reffe::GenericRefFE{RaviartThomas}, conf::DivConformity)
  get_face_dofs(reffe)
end

function get_dof_basis(reffe::GenericRefFE{RaviartThomas},phi::Field)
  cache = return_cache(get_dof_basis,reffe,phi)
  evaluate!(cache,get_dof_basis,reffe,phi)
end

function return_cache(::typeof(get_dof_basis),reffe::GenericRefFE{RaviartThomas},phi::Field)
  p = get_polytope(reffe)
  prebasis = get_prebasis(reffe)
  order = get_order(prebasis)
  et = return_type(prebasis)
  nf_nodes, nf_moments = _RT_nodes_and_moments(et,p,order,GenericField(identity))
  db = MomentBasedDofBasis(nf_nodes, nf_moments)

  face_moments = deepcopy(nf_moments)
  Jt_q_cache = return_cache(∇(phi),db.nodes)
  cache = (db.nodes, db.face_nodes, nf_moments, face_moments, Jt_q_cache)
  cache
end


function evaluate!(cache,::typeof(get_dof_basis),reffe::GenericRefFE{RaviartThomas},phi::Field)
  nodes, nf_nodes, nf_moments, face_moments, Jt_q_cache = cache

  Jt_q = evaluate!(Jt_q_cache,∇(phi),nodes)
  for face in 1:length(face_moments)
    moments = nf_moments[face]
    if length(moments) > 0
      num_qpoints, num_moments = size(moments)
      for i in 1:num_qpoints
        Jt_q_i = Jt_q[nf_nodes[face][i]]
        change = det(Jt_q_i) * inv(transpose(Jt_q_i))
        for j in 1:num_moments
          face_moments[face][i,j] = nf_moments[face][i,j] ⋅ change
        end
      end
    end
  end
  MomentBasedDofBasis(nodes,face_moments,nf_nodes)
end

function _RT_nodes_and_moments(::Type{et}, p::Polytope, order::Integer, phi::Field) where et

  D = num_dims(p)
  ft = VectorValue{D,et}
  pt = Point{D,et}

  nf_nodes = [ zeros(pt,0) for face in 1:num_faces(p)]
  nf_moments = [ zeros(ft,0,0) for face in 1:num_faces(p)]

  fcips, fmoments = _RT_face_values(p,et,order,phi)
  frange = get_dimrange(p,D-1)
  nf_nodes[frange] = fcips
  nf_moments[frange] = fmoments

  if (order > 0)
    ccips, cmoments = _RT_cell_values(p,et,order,phi)
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
  fgeomap = lazy_map(linear_combination,cfvs,cfshfs)
end

function _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)
  nc = length(fgeomap)
  c_fips = fill(fips,nc)
  c_wips = fill(wips,nc)
  pquad = lazy_map(evaluate,fgeomap,c_fips)

  if is_n_cube(p)
    c_detwips = c_wips
  elseif is_simplex(p)
    # Must account for diagonals in simplex discretizations to get the correct
    # scaling
    Jt1 = lazy_map(∇,fgeomap)
    Jt1_ips = lazy_map(evaluate,Jt1,c_fips)
    det_J = lazy_map(Broadcasting(meas),Jt1_ips)

    c_detwips = collect(lazy_map(Broadcasting(*),c_wips,det_J))
  end

  c_fips, pquad, c_detwips
end

function _broadcast(::Type{T},n,b) where T
  c = Array{T}(undef,size(b))
  for (ii, i) in enumerate(b)
    c[ii] = i⋅n
  end
  return c
end

function _RT_face_moments(p, fshfs, c_fips, fcips, fwips,phi)
  nc = length(c_fips)
  cfshfs = fill(fshfs, nc)
  cvals = lazy_map(evaluate,cfshfs,c_fips)
  cvals = [fwips[i].*cvals[i] for i in 1:nc]
  # fns, os = get_facet_normal(p)
  fns = get_facet_normal(p)
  os = get_facet_orientations(p)

  # Must express the normal in terms of the real/reference system of
  # coordinates (depending if phi≡I or phi is a mapping, resp.)
  # Hence, J = transpose(grad(phi))

  Jt = fill(∇(phi),nc)
  Jt_inv = lazy_map(Operation(inv),Jt)
  det_Jt = lazy_map(Operation(det),Jt)
  change = lazy_map(*,det_Jt,Jt_inv)
  change_ips = lazy_map(evaluate,change,fcips)

  cvals = [ _broadcast(typeof(n),n*o,J.*b) for (n,o,b,J) in zip(fns,os,cvals,change_ips)]

  return cvals
end

# It provides for every face the nodes and the moments arrays
function _RT_face_values(p,et,order,phi)

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
  fmoments = _RT_face_moments(p, fshfs, c_fips, fcips, fwips, phi)

  return fcips, fmoments

end

function _RT_cell_moments(p, cbasis, ccips, cwips)
  # Interior DOFs-related basis evaluated at interior integration points
  ishfs_iips = evaluate(cbasis,ccips)
  return cwips.⋅ishfs_iips
end

_p_filter(e,order) = (sum(e) <= order)

# It provides for every cell the nodes and the moments arrays
function _RT_cell_values(p,et,order,phi)
  # Compute integration points at interior
  degree = 2*(order+1)
  iquad = Quadrature(p,degree)
  ccips = get_coordinates(iquad)
  cwips = get_weights(iquad)

  # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
  if is_n_cube(p)
    cbasis = QGradMonomialBasis{num_dims(p)}(et,order-1)
  elseif is_simplex(p)
    T = VectorValue{num_dims(p),et}
    cbasis = MonomialBasis{num_dims(p)}(T,order-1, _p_filter)
  else
    @notimplemented
  end
  cell_moments = _RT_cell_moments(p, cbasis, ccips, cwips )

  # Must scale weights using phi map to get the correct integrals
  # scaling = det(grad(phi))
  Jt = ∇(phi)
  Jt_inv = inv(Jt)
  det_Jt = det(Jt)
  change = det_Jt*Jt_inv
  change_ips = evaluate(change,ccips)

  cmoments = change_ips.⋅cell_moments

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

function evaluate!(
  cache,
  ::Broadcasting{typeof(∇)},
  a::Fields.BroadcastOpFieldArray{ContraVariantPiolaMap})
  v, Jt, detJ = a.args
  # Assuming J comes from an affine map
  ∇v = Broadcasting(∇)(v)
  k = ContraVariantPiolaMap()
  Broadcasting(Operation(k))(∇v,Jt,detJ)
end

function lazy_map(
  ::Broadcasting{typeof(gradient)},
  a::LazyArray{<:Fill{Broadcasting{Operation{ContraVariantPiolaMap}}}})
  v, Jt, detJ = a.args
  ∇v = lazy_map(Broadcasting(∇),v)
  k = ContraVariantPiolaMap()
  lazy_map(Broadcasting(Operation(k)),∇v,Jt,detJ)
end

function evaluate!(cache,::ContraVariantPiolaMap,v::Number,Jt::Number,detJ::Number)
  v⋅((1/detJ)*Jt)
end

function evaluate!(cache,k::ContraVariantPiolaMap,v::AbstractVector{<:Field},phi::Field)
  Jt = ∇(phi)
  detJ = Operation(det)(Jt)
  Broadcasting(Operation(k))(v,Jt,detJ)
end

function lazy_map(
  k::ContraVariantPiolaMap,
  cell_ref_shapefuns::AbstractArray{<:AbstractArray{<:Field}},
  cell_map::AbstractArray{<:Field})

  cell_Jt = lazy_map(∇,cell_map)
  cell_detJ = lazy_map(Operation(det),cell_Jt)

  lazy_map(Broadcasting(Operation(k)),cell_ref_shapefuns,cell_Jt,cell_detJ)
end

PushForwardMap(reffe::GenericRefFE{RaviartThomas}) = ContraVariantPiolaMap()
