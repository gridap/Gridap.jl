struct BDM <: DivConforming end

const bdm = BDM()

"""
BDMRefFE(::Type{et},p::Polytope,order::Integer) where et

The `order` argument has the following meaning: the divergence of the  functions in this basis
is in the P space of degree `order-1`.

"""
function BDMRefFE(::Type{et},p::Polytope,order::Integer) where et

  D = num_dims(p)

  vet = VectorValue{num_dims(p),et}

  if is_simplex(p)
    prebasis = MonomialBasis(vet,p,order)
  else
    @notimplemented "BDM Reference FE only available for simplices"
  end

  nf_nodes, nf_moments = _BDM_nodes_and_moments(et,p,order,GenericField(identity))

  face_own_dofs = _face_own_dofs_from_moments(nf_moments)

  face_dofs = face_own_dofs

  dof_basis = MomentBasedDofBasis(nf_nodes, nf_moments)

  ndofs = num_dofs(dof_basis)

  metadata = nothing

  reffe = GenericRefFE{BDM}(
  ndofs,
  p,
  prebasis,
  dof_basis,
  DivConformity(),
  metadata,
  face_dofs)

  reffe
end

function ReferenceFE(p::Polytope,::BDM, order)
  BDMRefFE(Float64,p,order)
end

function ReferenceFE(p::Polytope,::BDM,::Type{T}, order) where T
  BDMRefFE(T,p,order)
end

function Conformity(reffe::GenericRefFE{BDM},sym::Symbol)
  hdiv = (:Hdiv,:HDiv)
  if sym == :L2
    L2Conformity()
  elseif sym in hdiv
    DivConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a BDM reference FE.

    Possible values of conformity for this reference fe are $((:L2, hdiv...)).
      """
    end
  end

  function get_face_own_dofs(reffe::GenericRefFE{BDM}, conf::DivConformity)
    get_face_dofs(reffe)
  end

  function _BDM_nodes_and_moments(::Type{et}, p::Polytope, order::Integer, phi::Field) where et

    D = num_dims(p)
    ft = VectorValue{D,et}
    pt = Point{D,et}

    nf_nodes = [ zeros(pt,0) for face in 1:num_faces(p)]
    nf_moments = [ zeros(ft,0,0) for face in 1:num_faces(p)]

    fcips, fmoments = _BDM_face_values(p,et,order,phi)
    frange = get_dimrange(p,D-1)
    nf_nodes[frange] = fcips
    nf_moments[frange] = fmoments

    if (order > 1)
      ccips, cmoments = _BDM_cell_values(p,et,order,phi)
      crange = get_dimrange(p,D)
      nf_nodes[crange] = ccips
      nf_moments[crange] = cmoments
    end

    nf_nodes, nf_moments
  end

  function _BDM_face_moments(p, fshfs, c_fips, fcips, fwips,phi)
    nc = length(c_fips)
    cfshfs = fill(fshfs, nc)
    cvals = lazy_map(evaluate,cfshfs,c_fips)
    cvals = [fwips[i].*cvals[i] for i in 1:nc]
    fns = get_facet_normal(p)

    # Must express the normal in terms of the real/reference system of
    # coordinates (depending if phi≡I or phi is a mapping, resp.)
    # Hence, J = transpose(grad(phi))

    Jt = fill(∇(phi),nc)
    Jt_inv = lazy_map(Operation(pinvJt),Jt)
    det_Jt = lazy_map(Operation(meas),Jt)
    change = lazy_map(*,det_Jt,Jt_inv)
    change_ips = lazy_map(evaluate,change,fcips)

    cvals = [ _broadcast(typeof(n),n,J.*b) for (n,b,J) in zip(fns,cvals,change_ips)]

    return cvals
  end

  # It provides for every face the nodes and the moments arrays
  function _BDM_face_values(p,et,order,phi)

    # Reference facet
    @check is_simplex(p) "We are assuming that all n-faces of the same n-dim are the same."
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
    degree = (order)*2
    fquad = Quadrature(fp,degree)
    fips = get_coordinates(fquad)
    wips = get_weights(fquad)

    c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

    # Moments (fmoments)
    # The BDM prebasis is expressed in terms of shape function
    fshfs = MonomialBasis(et,fp,order)

    # Face moments, i.e., M(Fi)_{ab} = q_RF^a(xgp_RFi^b) w_Fi^b n_Fi ⋅ ()
    fmoments = _BDM_face_moments(p, fshfs, c_fips, fcips, fwips, phi)

    return fcips, fmoments

  end

  function _BDM_cell_moments(p, cbasis, ccips, cwips)
    # Interior DOFs-related basis evaluated at interior integration points
    ishfs_iips = evaluate(cbasis,ccips)
    return cwips.⋅ishfs_iips
  end

  # It provides for every cell the nodes and the moments arrays
  function _BDM_cell_values(p,et,order,phi)
    # Compute integration points at interior
    degree = 2*(order)
    iquad = Quadrature(p,degree)
    ccips = get_coordinates(iquad)
    cwips = get_weights(iquad)

    # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
    if is_simplex(p)
      T = VectorValue{num_dims(p),et}
      # cbasis = GradMonomialBasis{num_dims(p)}(T,order-1)
      cbasis = Polynomials.NedelecPrebasisOnSimplex{num_dims(p)}(order-2)
    else
      @notimplemented
    end
    cell_moments = _BDM_cell_moments(p, cbasis, ccips, cwips )

    # Must scale weights using phi map to get the correct integrals
    # scaling = meas(grad(phi))
    Jt = ∇(phi)
    Jt_inv = pinvJt(Jt)
    det_Jt = meas(Jt)
    change = det_Jt*Jt_inv
    change_ips = evaluate(change,ccips)

    cmoments = change_ips.⋅cell_moments

    return [ccips], [cmoments]

  end
