"""
It return the `GenericRefFE` type for a Nedelec reference FE with order `order`
on the reference polytope `p`, using as scalar value `et`, i.e., `et` can be
`Float64` or `Float32`.
"""
# @santiagobadia : Project, go to complex numbers
function NedelecRefFE(p:: Polytope, et, order::Int)

  if !(all(extrusion(p).array .== HEX_AXIS))
    @notimplemented
  end

  # 1. Prebasis
  prebasis = GradMonomialBasis(VectorValue{dim(p),et},order)

  # Nface nodes, moments, and prebasis evaluated at nodes
  nf_nodes, nf_moments, pb_moments = _initialize_arrays(prebasis,p)

  # Face nodes and moments
  fcips, fmoments = _Nedelec_edge_values(p,et,order)

  # Insert nodes and moments in the nface arrays, together with the
  # prebasis evaluated at this nodes, for all edges
  nf_nodes,nf_moments,pb_moments = _insert_nface_values!(nf_nodes,nf_moments,pb_moments,prebasis,fcips,fmoments,p,1)

  # Idem for face values
  if ( dim(p) == 3 && order > 1)

    fcips, fmoments = _Nedelec_face_values(p,et,order)
    nf_nodes,nf_moments,pb_moments = _insert_nface_values!(nf_nodes,nf_moments,pb_moments,prebasis,fcips,fmoments,p,dim(p)-1)

  end

  # Idem for cell values
  if (order > 1)

    ccips, cmoments = _Nedelec_cell_values(p,et,order)
    nf_nodes,nf_moments,pb_moments = _insert_nface_values!(nf_nodes,nf_moments,pb_moments,prebasis,ccips,cmoments,p,dim(p))

  end

  _GenericRefFE(p,prebasis,nf_nodes,nf_moments,pb_moments)

end

# It provides for every edge the nodes and the moments arrays
function _Nedelec_edge_values(p,et,order)

  # Reference facet
  dim1 = 1
  fp = ref_nface_polytope(p,dim1) # dim 1

  # geomap from ref face to polytope faces
  fgeomap = _ref_face_to_faces_geomap(p,fp)

  # Compute integration points at all polynomial edges
  degree = order*2
  fquad = Quadrature(fp,degree)
  fips = coordinates(fquad)
  wips = weights(fquad)

  c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

  # Edge moments, i.e., M(Ei)_{ab} = q_RE^a(xgp_REi^b) w_Fi^b t_Ei ⋅ ()
  fshfs = Gridap.RefFEs._monomial_basis(fp,et,order-1)
  fmoments = _Nedelec_edge_moments(p, fshfs, c_fips, fcips, fwips)

  return fcips, fmoments

end

function _Nedelec_edge_moments(p, fshfs, c_fips, fcips, fwips)
  ts = edge_tangents(p)
  nc = length(c_fips)
  cfshfs = ConstantCellValue(fshfs, nc)
  cvals = evaluate(cfshfs,c_fips)
  cvals = [fwips[i]'.*cvals[i] for i in 1:nc]
  # @santiagobadia : Only working for oriented meshes now
  cvals = [ _broadcast(typeof(t),t,b) for (t,b) in zip(ts,cvals)]
  return cvals
end

# It provides for every face the nodes and the moments arrays
function _Nedelec_face_values(p,et,order)

  # Reference facet
  dimf = dim(p)-1
  fp = ref_nface_polytope(p,dimf) # dim 1

  # geomap from ref face to polytope faces
  fgeomap = _ref_face_to_faces_geomap(p,fp)

  # Compute integration points at all polynomial edges
  degree = order*2
  fquad = Quadrature(fp,degree)
  fips = coordinates(fquad)
  wips = weights(fquad)

  c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

  # Face moments, i.e., M(Fi)_{ab} = w_Fi^b q_RF^a(xgp_RFi^b) (n_Fi × ())
  fshfs = GradMonomialBasis(VectorValue{dim(fp),et},order-1)
  fmoments = _Nedelec_face_moments(p, fshfs, c_fips, fcips, fwips)

  return fcips, fmoments

end

function _Nedelec_face_moments(p, fshfs, c_fips, fcips, fwips)
  nc = length(c_fips)
  fvs = nfaces_vertices(p,dim(p)-1)
  fts = [hcat([vs[2]-vs[1]...],[vs[3]-vs[1]...]) for vs in fvs]
  cfshfs = ConstantCellValue(fshfs, nc)
  # Ref facet FE functions evaluated at the facet integration points (in ref facet)
  cvals = evaluate(cfshfs,c_fips)
  cvals = [fwips[i]'.*cvals[i] for i in 1:nc]
  fns, os = face_normals(p)
  # @santiagobadia : Temporary hack for making it work for structured hex meshes
  ft = eltype(fns)
  cvals = [ Gridap.RefFEs.GenericRefFEs._broadcast_extend(ft,Tm,b) for (Tm,b) in zip(fts,cvals)]
  cvals = [ Gridap.RefFEs.GenericRefFEs._broadcast_cross(ft,n*o,b) for (n,o,b) in zip(fns,os,cvals)]
  return cvals
end

# It provides for every cell the nodes and the moments arrays
function _Nedelec_cell_values(p,et,order)
  # Compute integration points at interior
  degree = 2*order
  iquad = Quadrature(p,degree)
  ccips = coordinates(iquad)
  cwips = weights(iquad)

  # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
  cbasis = CurlGradMonomialBasis(VectorValue{dim(p),et},order-1)
  cmoments = _Nedelec_cell_moments(p, cbasis, ccips, cwips )

  return [ccips], [cmoments]

end

const _Nedelec_cell_moments = _RT_cell_moments
