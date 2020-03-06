
"""
    RaviartThomasRefFE(::Type{et},p::Polytope,order::Integer) where et
"""
function RaviartThomasRefFE(::Type{et},p::Polytope,order::Integer) where et

  D = num_dims(p)

  prebasis = QCurlGradMonomialBasis{D}(et,order)

  nf_nodes, nf_moments = _RT_nodes_and_moments(et,p,order)

  face_own_dofs = _face_own_dofs_from_moments(nf_moments)

  face_own_dofs_permutations = _trivial_face_own_dofs_permutations(face_own_dofs)

  face_dofs = face_own_dofs

  dof_basis = MomentBasedDofBasis(nf_nodes, nf_moments)

  ndofs = num_dofs(dof_basis)

  reffe = GenericRefFE(
    ndofs,
    p,
    prebasis,
    dof_basis,
    face_own_dofs,
    face_own_dofs_permutations,
    face_dofs)

  reffe
end

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

  if (order > 1)
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
  fgeomap = lincomb(cfshfs,cfvs)
end

function _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)
  nc = length(fgeomap)
  c_fips = fill(fips,nc)
  c_wips = fill(wips,nc)
  pquad = evaluate(fgeomap,c_fips)
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
  cvals = evaluate(cfshfs,c_fips)
  cvals = [fwips[i].*cvals[i] for i in 1:nc]
  # fns, os = get_facet_normals(p)
  fns = get_facet_normals(p)
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
  degree = order*2
  fquad = Quadrature(fp,degree)
  fips = get_coordinates(fquad)
  wips = get_weights(fquad)

  c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

  # Moments (fmoments)
  # The RT prebasis is expressed in terms of shape function
  fshfs = MonomialBasis(et,fp,order-1)

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
  degree = 2*order
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

function _trivial_face_own_dofs_permutations(face_own_dofs)
  [ [collect(Int,1:length(dofs)),]  for dofs in face_own_dofs  ]
end

struct MomentBasedDofBasis{P,V} <: Dof
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

putaspasa() = "Hola"
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

function dof_cache(b::MomentBasedDofBasis,field)
  cf = field_cache(field,b.nodes)
  vals = evaluate_field!(cf,field,b.nodes)
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

function evaluate_dof!(cache,b::MomentBasedDofBasis,field)
  c, cf = cache
  vals = evaluate_field!(cf,field,b.nodes)
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
          dofs[o] += moments[i,j]*vals[nodes[i]]
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
            dofs[o,a] += moments[i,j]*vals[nodes[i],a]
          end
        end
        o += 1
      end
    end
  end
end
