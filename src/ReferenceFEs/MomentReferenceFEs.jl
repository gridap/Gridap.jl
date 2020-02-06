dim(p) = num_dims(p)

nfaces_dim(p,d) = collect(get_dimranges(p)[d+1])

function _initialize_arrays(prebasis,p)

  ft = VectorValue{num_dims(p),get_value_type(prebasis)}
  et = eltype(ft)
  pt = Point{num_dims(p),Float64}

  # Arrays of moments (as its evaluation for all prebasis shape functions)
  # and evaluation points per n-face
  nf_moments = Vector{Array{ft}}(undef,num_faces(p))
  nf_nodes = Vector{Array{pt}}(undef,num_faces(p))
  nshfs = num_terms(prebasis)
  pb_moments = zeros(et,nshfs,0)

  # Initialize to zero arrays
  zero_moments = zeros(ft,0,0)
  zero_nodes = zeros(pt,0)

  for dim in 0:num_dims(p)
    for inf in nfaces_dim(p,dim)
      nf_moments[inf] = zero_moments
      nf_nodes[inf] = zero_nodes
    end
  end

  return nf_nodes, nf_moments, pb_moments

end

"""
# Returns the reference polytope for n-faces of a given dimension.
"""
function ref_nface_polytope(p,nf_dim)
  nfs = nfaces_dim(p,nf_dim)
  fps = Gridap.ReferenceFEs._nface_ref_polytopes(p)[nfs]
  @assert(all(get_extrusion.(fps) .== get_extrusion(fps[1])), "All n-faces must be of the same type")
  return fps[1]
end

# Ref FE to faces geomaps
function _ref_face_to_faces_geomap(p,fp)
  cfvs = Gridap.ReferenceFEs._nfaces_vertices(p,dim(fp))
  nc = length(cfvs)
  freffe = LagrangianRefFE(Float64,fp,1)
  fshfs = get_shapefuns(freffe)
  cfshfs = fill(fshfs, nc)
  fgeomap = lincomb(cfshfs,cfvs)
end

function _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)
  nc = length(fgeomap)
  c_fips = fill(fips,nc)
  c_wips = fill(wips,nc)
  pquad = evaluate(fgeomap,c_fips)
  c_fips, pquad, c_wips
end

function _is_hex(p)
  for i in 2:length(get_extrusion(p))
    if get_extrusion(p)[i] != 1 return false && break end
  end
  return true
end

function _monomial_basis(p::Polytope{D}, T, order) where D
  if (_is_hex(p))
    # orders = order*ones(Int,D)
    prebasis = MonomialBasis{D}(T,order)
  elseif (_is_tet(p))
    filter(e,order) = sum(e) <= order
    prebasis = MonomialBasis{D}(T, filter, order)
  else
    @notimplemented
  end
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
  cvals = evaluate(cfshfs,c_fips)
  cvals = [fwips[i].*cvals[i] for i in 1:nc]
  # fns, os = get_facet_normals(p)
  fns = get_facet_normals(p)
  # @santiagobadia : Temporary hack for making it work for structured hex meshes
  # cvals = [ _broadcast(typeof(n),n*o,b) for (n,o,b) in zip(fns,os,cvals)]
  cvals = [ _broadcast(typeof(n),n,b) for (n,b) in zip(fns,cvals)]
  return cvals
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
  fips = get_coordinates(fquad)
  wips = get_weights(fquad)

  c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

  # Moments (fmoments)
  # The RT prebasis is expressed in terms of shape function
  fshfs = _monomial_basis(fp,et,order-1)


  # Face moments, i.e., M(Fi)_{ab} = q_RF^a(xgp_RFi^b) w_Fi^b n_Fi ⋅ ()
  fmoments = _RT_face_moments(p, fshfs, c_fips, fcips, fwips)

  return fcips, fmoments

end

function _nfaces_array_dim!(p,dim,array,nf_vals)
  nfs = nfaces_dim(p,dim)
  array[nfs] = nf_vals
end

function _insert_nface_values!(nf_nodes, nf_moments, pb_moments,
                               prebasis, fcips, fmoments, p, dim)

  # Evaluate basis in faces points, i.e., S(Fi)_{ab} = ϕ^a(xgp_Fi^b)
  pbasis_fcips = [evaluate(prebasis,ps) for ps in fcips]

  # Face moments evaluated for basis, i.e., DF = [S(F1)*M(F1)^T, …, S(Fn)*M(Fn)^T]
  fms_preb = [bps'*ms for (bps,ms) in zip(pbasis_fcips,fmoments)]

  _nfaces_array_dim!(p,dim,nf_moments,fmoments)
  _nfaces_array_dim!(p,dim,nf_nodes,fcips)
  pb_moments = hcat(pb_moments,fms_preb...)

  return nf_nodes, nf_moments, pb_moments

end

function _RT_cell_moments(p, cbasis, ccips, cwips)
  # Interior DOFs-related basis evaluated at interior integration points
  ishfs_iips = evaluate(cbasis,ccips)
  return cwips.*ishfs_iips
end

# It provides for every cell the nodes and the moments arrays
function _RT_cell_values(p,et,order)
  # Compute integration points at interior
  degree = 2*order
  iquad = Quadrature(p,degree)
  ccips = get_coordinates(iquad)
  cwips = get_weights(iquad)

  # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
  cbasis = QGradMonomialBasis{dim(p)}(et,order-1)
  cmoments = _RT_cell_moments(p, cbasis, ccips, cwips )

  return [ccips], [cmoments]

end

_num_nfaces(p) = length(get_faces(p))

function _nfacedofs_basis(p,moments)
  ndofs = [size(moments[i],1) for i in 1:length(moments)]
  _nfacedofs = Vector{Vector{Int}}(undef,_num_nfaces(p))
  _nfacedofs[1] = Int[1:ndofs[1]...]
  k = ndofs[1]
  for i in 2:length(ndofs)
    _nfacedofs[i] = [k+1:k+ndofs[i]...]
    k += ndofs[i]
  end
  return _nfacedofs
end

struct GenericDofBasis{P,V} <: Dof
  nodes::Vector{P}
  face_moments::Vector{Array{V}}
  face_nodes
end

function _GenericDofBasis(f_nodes,f_moments)
  _face_nodes = [length(i) for i in f_nodes]
  prepend!(_face_nodes,1)
  for i in 2:length(_face_nodes)
    _face_nodes[i] += _face_nodes[i-1]
  end
  _face_nodes
  face_nodes = [_face_nodes[i]:_face_nodes[i+1]-1 for i in 1:length(_face_nodes)-1]
  nodes = vcat([i for i in f_nodes]...)
  GenericDofBasis(nodes,f_moments,face_nodes)
end

function dof_cache(b::GenericDofBasis,field)
  cf = field_cache(field,b.nodes)
  vals = evaluate_field!(cf,field,b.nodes)
  ndofs = length(b.nodes)
  r = _generic_dof_cache(vals,ndofs)
  c = CachedArray(r)
  (c, cf)
end

function _generic_dof_cache(vals::AbstractVector,ndofs)
  T = eltype(vals)
  r = zeros(eltype(T),ndofs)
end

function _generic_dof_cache(vals::AbstractMatrix,ndofs)
  _, npdofs = size(vals)
  T = eltype(vals)
  r = zeros(eltype(T),ndofs,npdofs)
end

function evaluate_dof!(cache,b::GenericDofBasis,field)
  c, cf = cache
  vals = evaluate_field!(cf,field,b.nodes)
  k = 0
  dofs = c.array
  for (i,m) in enumerate(b.face_moments)
    l = size(m,2)
    if length(m) > 0
      dofs[k+1:k+l] = m'*vals[b.face_nodes[i]]
      k += l
    end
    println(i,m)
  end
  return dofs
end
