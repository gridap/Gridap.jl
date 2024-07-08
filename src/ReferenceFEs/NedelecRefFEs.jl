struct CurlConformity <: Conformity end

struct Nedelec <: ReferenceFEName end

const nedelec = Nedelec()

"""
    NedelecRefFE(::Type{et},p::Polytope,order::Integer) where et

The `order` argument has the following meaning: the curl of the  functions in this basis
is in the Q space of degree `order`.
"""
function NedelecRefFE(::Type{et},p::Polytope,order::Integer) where et

  # @santiagobadia : Project, go to complex numbers
  D = num_dims(p)

  if is_n_cube(p)
    prebasis = QGradMonomialBasis{D}(et,order)
  elseif is_simplex(p)
    prebasis = Polynomials.NedelecPrebasisOnSimplex{D}(order)
  else
    @unreachable "Only implemented for n-cubes and simplices"
  end

  nf_nodes, nf_moments = _Nedelec_nodes_and_moments(et,p,order)

  face_own_dofs = _face_own_dofs_from_moments(nf_moments)

  face_dofs = face_own_dofs

  dof_basis = MomentBasedDofBasis(nf_nodes, nf_moments)

  ndofs = num_dofs(dof_basis)

  metadata = nothing

  reffe = GenericRefFE{Nedelec}(
    ndofs,
    p,
    prebasis,
    dof_basis,
    CurlConformity(),
    metadata,
    face_dofs)

  reffe
end

function ReferenceFE(p::Polytope,::Nedelec, order)
  NedelecRefFE(Float64,p,order)
end

function ReferenceFE(p::Polytope,::Nedelec,::Type{T}, order) where T
  NedelecRefFE(T,p,order)
end

function Conformity(reffe::GenericRefFE{Nedelec},sym::Symbol)
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

function get_face_own_dofs(reffe::GenericRefFE{Nedelec}, conf::CurlConformity)
  reffe.face_dofs # For Nedelec, this member variable holds the face owned dofs
end

function get_face_own_dofs(reffe::GenericRefFE{Nedelec}, conf::L2Conformity)
  face_own_dofs=[Int[] for i in 1:num_faces(reffe)]
  face_own_dofs[end]=collect(1:num_dofs(reffe))
  face_own_dofs
end

function get_face_dofs(reffe::GenericRefFE{Nedelec,Dc}) where Dc
  face_dofs=[Int[] for i in 1:num_faces(reffe)]
  face_own_dofs=get_face_own_dofs(reffe)
  p = get_polytope(reffe)
  for d=1:Dc # Starting from edges, vertices do not own DoFs for Nedelec
    first_face = get_offset(p,d)
    nfaces     = num_faces(reffe,d)
    for face=first_face+1:first_face+nfaces
      for df=1:d-1
        face_faces  = get_faces(p,d,df)
        first_cface = get_offset(p,df)
        for cface in face_faces[face-first_face]
          cface_own_dofs = face_own_dofs[first_cface+cface]
          for dof in cface_own_dofs
              push!(face_dofs[face],dof)
          end
        end 
      end
      for dof in face_own_dofs[face]
        push!(face_dofs[face],dof)
      end
    end
  end
  face_dofs
end


function _Nedelec_nodes_and_moments(::Type{et}, p::Polytope, order::Integer) where et

  @notimplementedif !( is_n_cube(p) || (is_simplex(p) ) )

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

    if is_n_cube(p)
      fcips, fmoments = _Nedelec_face_values(p,et,order)
    else
      fcips, fmoments = _Nedelec_face_values_simplex(p,et,order)
    end

    frange = get_dimrange(p,D-1)
    nf_nodes[frange] = fcips
    nf_moments[frange] = fmoments

  end

  if ( is_n_cube(p) && order > 0) || ( is_simplex(p) && order > D-2)

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
  degree = (order)*2
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
  @assert is_n_cube(p) "We are assuming that all n-faces of the same n-dim are the same."
  fp = Polytope{num_dims(p)-1}(p,1)

  # geomap from ref face to polytope faces
  fgeomap = _ref_face_to_faces_geomap(p,fp)

  # Compute integration points at all polynomial edges
  degree = (order)*2
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
  cvals = lazy_map(evaluate,cfshfs,c_fips)

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

function _Nedelec_face_values_simplex(p,et,order)

  # Reference facet
  @assert is_simplex(p) "We are assuming that all n-faces of the same n-dim are the same."
  fp = Polytope{num_dims(p)-1}(p,1)

  # geomap from ref face to polytope faces
  fgeomap = _ref_face_to_faces_geomap(p,fp)

  # Compute integration points at all polynomial edges
  degree = (order)*2
  fquad = Quadrature(fp,degree)
  fips = get_coordinates(fquad)
  wips = get_weights(fquad)

  c_fips, fcips, fwips, fJtips = _nfaces_evaluation_points_weights_with_jac(p, fgeomap, fips, wips)

  Df = num_dims(fp)
  fshfs = MonomialBasis{Df}(VectorValue{Df,et},order-1,(e,k)->sum(e)<=k)

  fmoments = _Nedelec_face_moments_simplex(p, fshfs, c_fips, fcips, fwips, fJtips)

  return fcips, fmoments

end

function _nfaces_evaluation_points_weights_with_jac(p, fgeomap, fips, wips)
  nc = length(fgeomap)
  c_fips = fill(fips,nc)
  c_wips = fill(wips,nc)
  pquad = lazy_map(evaluate,fgeomap,c_fips)
  ## Must account for diagonals in simplex discretizations to get the correct
  ## scaling
  Jt1 = lazy_map(∇,fgeomap)
  Jt1_ips = lazy_map(evaluate,Jt1,c_fips)
  #det_J = lazy_map(Broadcasting(meas),Jt1_ips)
  #c_detwips = collect(lazy_map(Broadcasting(*),c_wips,det_J))
  c_detwips = c_wips
  c_fips, pquad, c_detwips, Jt1_ips
end

function _Nedelec_face_moments_simplex(p, fshfs, c_fips, fcips, fwips, fJtips)
  nc = length(c_fips)
  cfshfs = fill(fshfs, nc)
  cfshfs_fips = lazy_map(evaluate,cfshfs,c_fips)
  function weigth(qij,Jti,wi)
    Ji = transpose(Jti)
    Ji⋅qij*wi
  end
  cvals = map(Broadcasting(weigth),cfshfs_fips,fJtips,fwips)
  return cvals
end

# It provides for every cell the nodes and the moments arrays
function _Nedelec_cell_values(p,et,order)

  # Compute integration points at interior
  degree = 2*(order)
  iquad = Quadrature(p,degree)
  ccips = get_coordinates(iquad)
  cwips = get_weights(iquad)

  # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
  if is_n_cube(p)
    cbasis = QCurlGradMonomialBasis{num_dims(p)}(et,order-1)
  else
    D = num_dims(p)
    cbasis = MonomialBasis{D}(VectorValue{D,et},order-D+1,(e,k)->sum(e)<=k)
  end
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

struct CoVariantPiolaMap <: Map end

function evaluate!(
  cache,
  ::Broadcasting{typeof(∇)},
  a::Fields.BroadcastOpFieldArray{CoVariantPiolaMap})
  v, Jt = a.args
  # Assuming J comes from an affine map
  ∇v = Broadcasting(∇)(v)
  k = CoVariantPiolaMap()
  Broadcasting(Operation(k))(∇v,Jt)
end

function lazy_map(
  ::Broadcasting{typeof(gradient)},
  a::LazyArray{<:Fill{Broadcasting{Operation{CoVariantPiolaMap}}}})
  v, Jt = a.args
  ∇v = lazy_map(Broadcasting(∇),v)
  k = CoVariantPiolaMap()
  lazy_map(Broadcasting(Operation(k)),∇v,Jt)
end

function evaluate!(cache,::CoVariantPiolaMap,v::Number,Jt::Number)
  v⋅transpose(inv(Jt)) # we multiply by the right side to compute the gradient correctly
end

function evaluate!(cache,k::CoVariantPiolaMap,v::AbstractVector{<:Field},phi::Field)
  Jt = ∇(phi)
  Broadcasting(Operation(k))(v,Jt)
end

function lazy_map(
  k::CoVariantPiolaMap,
  cell_ref_shapefuns::AbstractArray{<:AbstractArray{<:Field}},
  cell_map::AbstractArray{<:Field})

  cell_Jt = lazy_map(∇,cell_map)
  lazy_map(Broadcasting(Operation(k)),cell_ref_shapefuns,cell_Jt)
end
