"""
It return the `GenericRefFE` type for a Raviart-Thomas reference FE with order `order`
on the reference polytope `p`, using as scalar value `et`, i.e., `et` can be
`Float64` or `Float32`.
"""
function RTRefFE(p:: Polytope, et, order::Int)

  if !(all(extrusion(p).array .== HEX_AXIS))
    @notimplemented
  end

  # 1. Prebasis
  prebasis = CurlGradMonomialBasis(VectorValue{dim(p),et},order)

  # Nface nodes, moments, and prebasis evaluated at nodes
  nf_nodes, nf_moments, pb_moments = _initialize_arrays(prebasis,p)

  # Face values
  fcips, fmoments = _RT_face_values(p,et,order)
  nf_nodes,nf_moments,pb_moments = _insert_nface_values!(nf_nodes,nf_moments,pb_moments,prebasis,fcips,fmoments,p,dim(p)-1)

  # Cell values
  if (order > 1)

    ccips, cmoments = _RT_cell_values(p,et,order)
    nf_nodes,nf_moments,pb_moments = _insert_nface_values!(nf_nodes,nf_moments,pb_moments,prebasis,ccips,cmoments,p,dim(p))

  end

  _GenericRefFE(p,prebasis,nf_nodes,nf_moments,pb_moments)

end

# It provides for every face the nodes and the moments arrays
function _RT_face_values(p,et,order)

  # Reference facet
  fp = ref_nface_polytope(p,dim(p)-1)

  # geomap from ref face to polytope faces
  fgeomap = _ref_face_to_faces_geomap(p,fp)

  # Nodes are integration points (for exact integration)
  # Thus, we define the integration points in the reference
  # face polytope (fips and wips). Next, we consider the
  # n-face-wise arrays of nodes in fp (constant cell array c_fips)
  # the one of the points in the polytope after applying the geopmap
  # (fcips), and the weights for these nodes (fwips, a constant cell array)
  # Nodes (fcips)
  degree = order*2
  fquad = Quadrature(fp,degree)
  fips = coordinates(fquad)
  wips = weights(fquad)

  c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

  # Moments (fmoments)
  # The RT prebasis is expressed in terms of shape function
  fshfs = Gridap.RefFEs._monomial_basis(fp,et,order-1)

  # Face moments, i.e., M(Fi)_{ab} = q_RF^a(xgp_RFi^b) w_Fi^b n_Fi ⋅ ()
  fmoments = _RT_face_moments(p, fshfs, c_fips, fcips, fwips)

  return fcips, fmoments

end

function _RT_face_moments(p, fshfs, c_fips, fcips, fwips)
  nc = length(c_fips)
  cfshfs = ConstantCellValue(fshfs, nc)
  cvals = evaluate(cfshfs,c_fips)
  cvals = [fwips[i]'.*cvals[i] for i in 1:nc]
  fns, os = face_normals(p)
  # @santiagobadia : Temporary hack for making it work for structured hex meshes
  cvals = [ _broadcast(typeof(n),n*o,b) for (n,o,b) in zip(fns,os,cvals)]
  return cvals
end

# It provides for every cell the nodes and the moments arrays
function _RT_cell_values(p,et,order)
    # Compute integration points at interior
    degree = 2*order
    iquad = Quadrature(p,degree)
    ccips = coordinates(iquad)
    cwips = weights(iquad)

    # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
    cbasis = GradMonomialBasis(VectorValue{dim(p),et},order-1)
    cmoments = _RT_cell_moments(p, cbasis, ccips, cwips )

    return [ccips], [cmoments]

  end

function _RT_cell_moments(p, cbasis, ccips, cwips)
  # Interior DOFs-related basis evaluated at interior integration points
  ishfs_iips = evaluate(cbasis,ccips)
  return cwips'.*ishfs_iips
end
