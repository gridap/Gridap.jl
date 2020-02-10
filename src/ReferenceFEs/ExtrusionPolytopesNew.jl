
"""
Constant to be used in order to indicate a hex-like extrusion axis.
"""
const HEX_AXIS = 1

"""
Constant to be used in order to indicate a tet-like extrusion axis.
"""
const TET_AXIS = 2

struct NFace{D}
  anchor::Point{D,Int}
  extrusion::Point{D,Int}
end

struct DFace{D} <: GridapType
  extrusion::Point{D,Int}
  nfaces::Vector{NFace{D}}
  dimranges::Vector{UnitRange{Int}}
  dims::Vector{Int}
  nf_nfs::Vector{Vector{Int}}
  nf_dimranges::Vector{Vector{UnitRange{Int64}}}
  nf_dims::Vector{Vector{Int}}
end

function DFace(extrusion::Point{D,Int}) where D
  zerop = zero(Point{D,Int})
  pol_nfs, pol_dim = _polytopenfaces(zerop, extrusion)
  nfs_id = Dict{NFace{D},Int}(nf => i for (i, nf) in enumerate(pol_nfs))
  nf_nfs, nf_dimranges = _polytopemesh(pol_nfs, nfs_id)
  dims = _nface_to_nfacedim(length(pol_nfs), pol_dim)
  nf_dims = _nf_dims(nf_nfs,nf_dimranges)
  DFace{D}(extrusion, pol_nfs, pol_dim, dims, nf_nfs, nf_dimranges, nf_dims)
end

# Helpers for the construction of a DFace

# Provides for all n-faces of a polytope the d-faces for 0 <= d <n on its
# boundary (e.g., given a face, it provides the set of edges and corners on its
# boundary) using the global n-face numbering of the base polytope
function _polytopemesh(nfaces::Vector{NFace{D}}, nfaceid::Dict{NFace{D},Int}) where D
  num_nfs = length(nfaces)
  nfnfs = Vector{Vector{Int}}(undef, num_nfs)
  nfnfs_dim = Vector{Vector{UnitRange{Int}}}(undef, num_nfs)
  for (inf, nf) in enumerate(nfaces)
    nf_nfs, dimnfs = _polytopenfaces(nf.anchor, nf.extrusion)
    nfnfs[inf] = [ nfaceid[nf] for nf in nf_nfs]
    nfnfs_dim[inf] = dimnfs
  end

  return (nfnfs, nfnfs_dim)
end

function _nface_to_nfacedim(nnfaces::Int,dimranges::Vector{UnitRange{Int}})
  v = zeros(Int,nnfaces)
  for (d,r) in enumerate(dimranges)
    v[r] .= d-1
  end
  v
end

function _nf_dims(nf_nfs,nf_dimranges)
  nf_dims = Vector{Int}[]
  for i in 1:length(nf_nfs)
    nfs = nf_nfs[i]
    dimranges = nf_dimranges[i]
    dims = _nface_to_nfacedim(length(nfs),dimranges)
    push!(nf_dims,dims)
  end
  nf_dims
end

# Generates the array of n-faces of a polytope
function _polytopenfaces(anchor, extrusion)

  D = n_components(extrusion)
  zerop = zero(Point{D,Int})
  nf_nfs = Vector{NFace{D}}(undef,0)
  _nfaceboundary!(anchor, zerop, extrusion, true, nf_nfs)

  for i in 1:D
    sort!(nf_nfs, by = x -> x.anchor[i])
  end
  for i in 1:D
    sort!(nf_nfs, by = x -> x.extrusion[i])
  end
  for i in 1:D
    sort!(nf_nfs, by = x -> sum(x.extrusion))
  end

  numnfs = length(nf_nfs)
  nfsdim = [_nfdim(nf_nfs[i].extrusion) for i = 1:numnfs]
  dnf = _nfdim(extrusion)
  dimnfs = Array{UnitRange{Int64},1}(undef, dnf + 1)
  dim = 0
  i = 1
  for iface = 1:numnfs
    if (nfsdim[iface] > dim)
      # global dim; # global i
      dim += 1
      dimnfs[dim] = i:iface-1
      i = iface
    end
  end
  dimnfs[dnf+1] = numnfs:numnfs

  (nf_nfs, dimnfs)
end

function _nfdim(extrusion)
  z = zero(eltype(extrusion))
  for e in extrusion
    if e > 0
      z +=1
    end
  end
  z
end

# Generates the list of n-face of a polytope the d-faces for 0 <= d <n on its
# boundary
function _nfaceboundary!(anchor, extrusion, extend, isanchor, list)

  D = n_components(extrusion)
  newext = extend
  push!(list,NFace(anchor,extrusion))

  for i = 1:D
    curex = newext[i]
    if (curex > 0) # Perform extension

      newext = _newext(newext,i)
      edim = _edim(newext,i)
      tetp = _tetp(anchor,i)

      if (curex == HEX_AXIS) # Quad extension
        list = _nfaceboundary!(anchor + edim, extrusion, newext, false, list)
      elseif isanchor
        list = _nfaceboundary!(tetp, extrusion, newext, false, list)
      end
      list = _nfaceboundary!( anchor, extrusion + edim * curex, newext, false, list)

    end
  end

  list
end

function _newext(newext,i)
  m = zero(mutable(newext))
  D = n_components(newext)
  for j in 1:D
    m[j] = j == i ? 0 : newext[j]
  end
  Point(m)
end

function _edim(newext,i)
  m = zero(mutable(newext))
  D = n_components(newext)
  for j in 1:D
    m[j] = j == i ? 1 : 0
  end
  Point(m)
end

function _tetp(anchor,i)
  m = zero(mutable(anchor))
  D = n_components(anchor)
  for j in 1:D
    m[j] = j >= i ? anchor[j] : 0
  end
  Point(m)
end

# Helper API on a DFace

function _dimfrom_fs_dimto_fs(p::DFace{D}, dim_from::Int, dim_to::Int) where D
  @assert dim_from <= D
  if dim_to == dim_from
    return [ [i,] for i in 1:_num_nfaces(p,dim_from)]
  end
  @assert dim_to < dim_from
  dim_from += 1
  dim_to += 1
  dffs_r = p.nf_dimranges[end][dim_from]
  dffs_dtfs = Vector{Vector{Int}}(undef, dffs_r[end] - dffs_r[1] + 1)
  offs = p.nf_dimranges[end][dim_to][1] - 1
  for (i_dff, dff) in enumerate(dffs_r)
    dff_nfs = p.nf_nfs[dff]
    dff_dtfs_r = p.nf_dimranges[dff][dim_to]
    dff_dtfs = dff_nfs[dff_dtfs_r]
    dffs_dtfs[i_dff] = dff_dtfs .- offs
  end
  dffs_dtfs
end

# iface in global numeration 
function DFace{d}(p::DFace{D},iface::Int) where {d,D}
  @assert d == p.dims[iface]
  nf = p.nfaces[iface]
  r_ext = _eliminate_zeros(Val{d}(),nf.extrusion)
  DFace(r_ext)
end

function DFace{D}(p::DFace{D},iface::Int) where D
  @assert D == p.dims[iface]
  p
end

function _eliminate_zeros(::Val{d},a) where d
  b = zero(mutable(Point{d,Int}))
  D = n_components(a)
  k = 1
  for i in 1:D
    m = a[i]
    if (m != 0)
      b[k] = m
      k += 1
    end
  end
  Point(b)
end

# It generates the list of coordinates of all vertices in the polytope. It is
# assumed that the polytope has the bounding box [0,1]**dim
function _vertices_coordinates(::Type{T},p::DFace{D}) where {D,T}
  vs = first(_dimfrom_fs_dimto_fs(p, D, 0))
  vcs = zeros(Point{D,T},length(vs))
  for i = 1:length(vs)
    vx = p.nfaces[vs[i]]
    vc = vx.anchor.array.data
    vcs[i] = vc
  end
  vcs
end

# Returns number of nfaces of dimension dim
function _num_nfaces(polytope::DFace, dim::Integer)
  n = length(polytope.nf_dimranges)
  k = 0
  for nface = 1:n
    d = length(polytope.nf_dimranges[nface]) - 1
    if d == dim
      k += 1
    end
  end
  k
end

# Return the n-faces vertices coordinates array for a given n-face dimension
function _nfaces_vertices(::Type{T},p::DFace,d::Integer) where T
  nc = _num_nfaces(p,d)
  verts = _vertices_coordinates(T,p)
  faces_vs = _dimfrom_fs_dimto_fs(p,d,0)
  cfvs = collect(LocalToGlobalArray(faces_vs,verts))
end

# It generates the outwards normals of the facets of a polytope. It returns two
# arrays, the first one being the outward normal and the second one the orientation.
function _face_normals(::Type{T},p::DFace{D}) where {D,T}
  nf_vs = _dimfrom_fs_dimto_fs(p, D - 1, 0)
  vs = _vertices_coordinates(T,p)
  n = length(p.nf_dimranges[end][end-1])
  f_ns = zeros(Point{D,T},n)
  f_os = zeros(Int,n)
  for i_f = 1:n
    f_n, f_o = _facet_normal(T,p, nf_vs, vs, i_f)
    f_ns[i_f] = f_n
    f_os[i_f] = f_o
  end
  f_ns, f_os
end

function _face_normals(::Type{T},p::DFace{0}) where T
  (VectorValue{0,T}[], Int[])
end

function _facet_normal(::Type{T},p::DFace{D}, nf_vs, vs, i_f) where {D,T}
  if length(p.extrusion) > 1
    n1 = length(nf_vs[i_f])-1
    n2 = D
    v = zeros(T,n1,n2)
    for i in 2:length(nf_vs[i_f])
      vi = vs[nf_vs[i_f][i]] - vs[nf_vs[i_f][1]]
      for d in 1:n_components(vi)
        v[i-1,d] = vi[d]
      end
    end
    n = nullspace(v)
    n = n .* 1 / sqrt(dot(n, n))
    ext_v = _vertex_not_in_facet(p, i_f, nf_vs)
    v3 = vs[nf_vs[i_f][1]] - vs[ext_v]
    f_or = 1
    if dot(v3, n) < 0.0
      n *= -1
      f_or = -1
    end
  elseif (length(p.extrusion) == 1)
    ext_v = _vertex_not_in_facet(p, i_f, nf_vs)
    _n = vs[nf_vs[i_f][1]] - vs[ext_v]
    _n = _n .* 1 / dot(_n, _n)
    n = _n.array
    f_or = 1
  else
    @unreachable "O-dim polytopes do not have properly define outward facet normals"
  end
  Point{D,T}(n), f_or
end

function _vertex_not_in_facet(p::DFace, i_f, nf_vs)
  for i in p.nf_dimranges[end][1]
    is_in_f = false
    for j in nf_vs[i_f]
      if i == j
        is_in_f = true
        break
      end
    end
    if !is_in_f
      return i
      break
    end
  end
  -1
end

# It generates the tangent vectors for polytope edges.
function _edge_tangents(::Type{T},p::DFace{D}) where {D,T}
  ed_vs = _nfaces_vertices(T,p,1)
  ts = [(t = vs[2]-vs[1])/norm(t) for vs in ed_vs ]
  ts
end

function _edge_tangents(::Type{T},p::DFace{0}) where T
  VectorValue{0,T}[]
end

# It generates all the admissible permutations of nodes that lead to an
# admissible polytope
function _admissible_permutations(p::DFace{D}) where D
  if D > 3
    @warn "Computing permutations for a polytope of dim > 3 is overkill"
  end
  if D==0 || all( p.extrusion.array.data .== TET_AXIS )
    perms = _admissible_permutations_simplex(p)
  elseif all( p.extrusion.array.data .== HEX_AXIS)
    perms = _admissible_permutations_n_cube(p)
  else
    @notimplemented "admissible vertex permutations only implemented for simplices and n-cubes"
  end
  perms
end

function _admissible_permutations_simplex(p::DFace{D}) where D
  vs = first(_dimfrom_fs_dimto_fs(p, D, 0))
  num_vs = length(vs)
  l = [i for i = 1:num_vs]
  perms = Combinatorics.permutations(l)
  collect(perms)
end

function _admissible_permutations_n_cube(p::DFace{D}) where D
  vs = first(_dimfrom_fs_dimto_fs(p, D, 0))
  num_vs = length(vs)
  l = [i for i = 1:num_vs]
  perms = Combinatorics.permutations(l)
  vertices = p.nfaces[first(p.dimranges)]
  permuted_vertices = similar(vertices)
  admissible_perms = Vector{Int}[]
  grads = _setup_aux_grads(vertices)
  vol = -1.0
  for perm in perms
    for (j,cj) in enumerate(perm)
      permuted_vertices[j] = vertices[cj]
    end
    jac = _setup_aux_jacobian(grads,permuted_vertices)
    vol_i = abs(det(jac))
    if vol < 0
      vol = vol_i
    end
    if vol_i == vol
      push!(admissible_perms,perm)
    end
  end
  admissible_perms
end

function _setup_aux_grads(vertices::Vector{NFace{D}}) where D
  grads = zeros(Point{D,Int},length(vertices))
  m = zero(mutable(Point{D,Int}))
  for (i,vertex) in enumerate(vertices)
    x = vertex.anchor
    for di in 1:D
      v = 1
      for k in 1:D
        if k == di && x[k] == 0
          v *= -1
        end
      end
      m[di] = v
    end
    grads[i] = m
  end
  grads
end

function _setup_aux_jacobian(grads,permuted_vertices::Vector{NFace{D}}) where D
  p0 = zero(Point{D,Int})
  m = zero(mutable(outer(p0,p0)))
  for (i,pvertex) in enumerate(permuted_vertices)
    x = pvertex.anchor
    g = grads[i]
    for di in 1:D
      v = g[di]
      for dj in 1:D
        m[di,dj] += x[dj]*v
      end
    end
  end
  TensorValue(m)
end


