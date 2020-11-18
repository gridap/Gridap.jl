struct CurlConformity <: Conformity end

"""
    NedelecRefFE(::Type{et},p::Polytope,order::Integer) where et

The `order` argument has the following meaning: the curl of the  functions in this basis
is in the Q space of degree `order`.
"""
function NedelecRefFE(::Type{et},p::Polytope,order::Integer) where et

  # @santiagobadia : Project, go to complex numbers
  D = num_dims(p)

  prebasis = QGradMonomialBasis{D}(et,order)

  nf_nodes, nf_moments = _Nedelec_nodes_and_moments(et,p,order)

  face_own_dofs = _face_own_dofs_from_moments(nf_moments)

  face_dofs = face_own_dofs

  dof_basis = MomentBasedDofBasis(nf_nodes, nf_moments)

  ndofs = num_dofs(dof_basis)

  metadata = nothing

  reffe = GenericRefFE{:Nedelec}(
    ndofs,
    p,
    prebasis,
    dof_basis,
    CurlConformity(),
    metadata,
    face_dofs)

  reffe
end

function ReferenceFE(p::Polytope,::Val{:Nedelec}, order)
  NedelecRefFE(Float64,p,order)
end

function ReferenceFE(p::Polytope,::Val{:Nedelec},::Type{T}, order) where T
  NedelecRefFE(T,p,order)
end

function Conformity(reffe::GenericRefFE{:Nedelec},sym::Symbol)
  hcurl = (:Hcurl,:HCurl)
  if sym == :L2
    L2Conformity()
  elseif sym in hcurl
    CurlConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a Nedelec reference FE.

    Possible values of conformity for this reference fe are $((:L2, hcurl...)).
    """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{:Nedelec}, conf::CurlConformity)
  get_face_dofs(reffe)
end

function _Nedelec_nodes_and_moments(::Type{et}, p::Polytope, order::Integer) where et

  @notimplementedif ! is_n_cube(p)

  D = num_dims(p)
  ft = VectorValue{D,et}
  pt = Point{D,et}

  nf_nodes = [ zeros(pt,0) for face in 1:num_faces(p)]
  nf_moments = [ zeros(ft,0,0) for face in 1:num_faces(p)]

  ecips, emoments = _Nedelec_edge_values(p,et,order)
  erange = get_dimrange(p,1)
  nf_nodes[erange] = ecips
  nf_moments[erange] = emoments

  if ( num_dims(p) == 3 && order > 0)

    fcips, fmoments = _Nedelec_face_values(p,et,order)
    frange = get_dimrange(p,D-1)
    nf_nodes[frange] = fcips
    nf_moments[frange] = fmoments

  end

  if (order > 0)

    ccips, cmoments = _Nedelec_cell_values(p,et,order)
    crange = get_dimrange(p,D)
    nf_nodes[crange] = ccips
    nf_moments[crange] = cmoments

  end

  nf_nodes, nf_moments
end

function _Nedelec_edge_values(p,et,order)

  # Reference facet
  dim1 = 1
  ep = Polytope{dim1}(p,1)

  # geomap from ref face to polytope faces
  egeomap = _ref_face_to_faces_geomap(p,ep)

  # Compute integration points at all polynomial edges
  degree = (order+1)*2
  equad = Quadrature(ep,degree)
  cips = get_coordinates(equad)
  wips = get_weights(equad)


  c_eips, ecips, ewips = _nfaces_evaluation_points_weights(p, egeomap, cips, wips)

  # Edge moments, i.e., M(Ei)_{ab} = q_RE^a(xgp_REi^b) w_Fi^b t_Ei ⋅ ()
  eshfs = MonomialBasis(et,ep,order)
  emoments = _Nedelec_edge_moments(p, eshfs, c_eips, ecips, ewips)

  return ecips, emoments

end

function _Nedelec_edge_moments(p, fshfs, c_fips, fcips, fwips)
  ts = get_edge_tangent(p)
  nc = length(c_fips)
  cfshfs = fill(fshfs, nc)
  cvals = lazy_map(evaluate,cfshfs,c_fips)
  cvals = [fwips[i].*cvals[i] for i in 1:nc]
  # @santiagobadia : Only working for oriented meshes now
  cvals = [ _broadcast(typeof(t),t,b) for (t,b) in zip(ts,cvals)]
  return cvals
end

function _Nedelec_face_values(p,et,order)

  # Reference facet
  @assert is_simplex(p) || is_n_cube(p) "We are assuming that all n-faces of the same n-dim are the same."
  fp = Polytope{num_dims(p)-1}(p,1)

  # geomap from ref face to polytope faces
  fgeomap = _ref_face_to_faces_geomap(p,fp)

  # Compute integration points at all polynomial edges
  degree = (order+1)*2
  fquad = Quadrature(fp,degree)
  fips = get_coordinates(fquad)
  wips = get_weights(fquad)

  c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

  # Face moments, i.e., M(Fi)_{ab} = w_Fi^b q_RF^a(xgp_RFi^b) (n_Fi × ())
  fshfs = QGradMonomialBasis{num_dims(fp)}(et,order-1)

  fmoments = _Nedelec_face_moments(p, fshfs, c_fips, fcips, fwips)

  return fcips, fmoments

end

function _Nedelec_face_moments(p, fshfs, c_fips, fcips, fwips)
  nc = length(c_fips)
  cfshfs = fill(fshfs, nc)
  cvals = evaluate(cfshfs,c_fips)

  fvs = _nfaces_vertices(Float64,p,num_dims(p)-1)
  fts = [hcat([vs[2]-vs[1]...],[vs[3]-vs[1]...]) for vs in fvs]

  # Ref facet FE functions evaluated at the facet integration points (in ref facet)
  cvals = [fwips[i].*cvals[i] for i in 1:nc]

  fns = get_facet_normal(p)
  os = get_facet_orientations(p)
  # @santiagobadia : Temporary hack for making it work for structured hex meshes
  ft = eltype(fns)
  cvals = [ _broadcast_extend(ft,Tm,b) for (Tm,b) in zip(fts,cvals)]
  cvals = [ _broadcast_cross(ft,n*o,b) for (n,o,b) in zip(fns,os,cvals)]
  return cvals
end

# It provides for every cell the nodes and the moments arrays
function _Nedelec_cell_values(p,et,order)

  # Compute integration points at interior
  degree = 2*(order+1)
  iquad = Quadrature(p,degree)
  ccips = get_coordinates(iquad)
  cwips = get_weights(iquad)

  # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
  cbasis = QCurlGradMonomialBasis{num_dims(p)}(et,order-1)
  cmoments = _Nedelec_cell_moments(p, cbasis, ccips, cwips )

  return [ccips], [cmoments]

end

const _Nedelec_cell_moments = _RT_cell_moments

function _broadcast_cross(::Type{T},n,b) where T
  c = Array{T}(undef,size(b))
  for (ii, i) in enumerate(b)
    c[ii] = T(cross(get_array(i),get_array(n)))# cross product
  end
  return c
end

function _broadcast_extend(::Type{T},Tm,b) where T
  c = Array{T}(undef,size(b))
  for (ii,i) in enumerate(b)
    c[ii] = T(Tm*[i...])
  end
  return c
end
