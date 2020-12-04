
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
  nf_dimranges::Vector{Vector{UnitRange{Int}}}
  nf_dims::Vector{Vector{Int}}
end

"""
    struct ExtrusionPolytope{D} <: Polytope{D}
      extrusion::Point{D,Int}
      # + private fields
    end

Concrete type for polytopes that can be
represented with an "extrusion" tuple. The underlying extrusion is available
in the field `extrusion`. Instances of this type can be obtained with the constructors

- [`Polytope(extrusion::Int...)`](@ref)
- [`ExtrusionPolytope(extrusion::Int...)`](@ref)

"""
struct ExtrusionPolytope{D} <: Polytope{D}
  extrusion::Point{D,Int}
  dface::DFace{D}
  vertex_coords::Vector{Point{D,Float64}}
  face_normals::Vector{Point{D,Float64}}
  face_orientations::Vector{Int}
  vertex_perms::Vector{Vector{Int}}
  face_vertex_perms::Vector{Vector{Vector{Int}}}
  m_n_to_mface_to_nface::Matrix{Vector{Vector{Int}}}
  function ExtrusionPolytope(dface::DFace{D}) where D
    vertex_coords = _vertices_coordinates(Float64,dface)
    face_normals, face_orientations = _face_normals(Float64,dface)
    vertex_perms = _precompute_vertex_perms_if_possible(dface)
    face_vertex_perms = _admissible_face_vertex_permutations(dface)
    m_n_to_mface_to_nface = _precompute_m_n_to_mface_to_nface(dface)
    new{D}(
      dface.extrusion,dface,vertex_coords,face_normals,
      face_orientations,vertex_perms,face_vertex_perms,
      m_n_to_mface_to_nface)
  end
end

# Constructors

"""
    Polytope(extrusion::Int...)

Equivalent to `ExtrusionPolytope(extrusion...)`
"""
Polytope(extrusion::Int...) = Polytope(extrusion)

function Polytope(extrusion::NTuple{D,Int}) where D
  ExtrusionPolytope(extrusion)
end

"""
    ExtrusionPolytope(extrusion::Int...)

Generates an `ExtrusionPolytope` from the tuple `extrusion`.
The values in `extrusion` are either equal to the constant
[`HEX_AXIS`](@ref) or the constant [`TET_AXIS`](@ref).

# Examples

Creating a quadrilateral, a triangle, and a wedge

```jldoctest
using Gridap.ReferenceFEs

quad = ExtrusionPolytope(HEX_AXIS,HEX_AXIS)

tri = ExtrusionPolytope(TET_AXIS,TET_AXIS)

wedge = ExtrusionPolytope(TET_AXIS,TET_AXIS,HEX_AXIS)

println(quad == QUAD)
println(tri == TRI)
println(wedge == WEDGE)

# output
true
true
true

```

"""
function ExtrusionPolytope(extrusion::Int...)
  ExtrusionPolytope(extrusion)
end

function ExtrusionPolytope(extrusion::NTuple{N,Int}) where N
  ExtrusionPolytope(Point{N,Int}(extrusion))
end

function ExtrusionPolytope(extrusion::Point{N,Int}) where N
  dface = DFace(extrusion)
  ExtrusionPolytope(dface)
end

# Implementation of the interface

get_faces(p::ExtrusionPolytope) = p.dface.nf_nfs

function get_faces(p::ExtrusionPolytope,dimfrom::Integer,dimto::Integer)
  p.m_n_to_mface_to_nface[dimfrom+1,dimto+1]
end

get_dimranges(p::ExtrusionPolytope) = p.dface.dimranges

get_facedims(p::ExtrusionPolytope) = p.dface.dims

function get_face_dimranges(p::ExtrusionPolytope)
  p.dface.nf_dimranges
end

function get_face_dimranges(p::ExtrusionPolytope,d::Integer)
  r = get_dimrange(p,d)
  p.dface.nf_dimranges[r]
end

function Polytope{D}(p::ExtrusionPolytope,Dfaceid::Integer) where D
  offset = get_offset(p,D)
  id = Dfaceid + offset
  dface = DFace{D}(p.dface,id)
  ExtrusionPolytope(dface)
end

function Polytope{0}(p::ExtrusionPolytope,Dfaceid::Integer)
  @assert Dfaceid in 1:num_vertices(p)
  VERTEX
end

function Polytope{D}(p::ExtrusionPolytope{D},Dfaceid::Integer) where D
  @assert Dfaceid == 1
  p
end

function (==)(a::ExtrusionPolytope{D},b::ExtrusionPolytope{D}) where D
  #The first axis is irrelevant here
  ea = Point(Tuple(a.extrusion)[2:end])
  eb = Point(Tuple(b.extrusion)[2:end])
  ea == eb
end

function (==)(a::ExtrusionPolytope{1},b::ExtrusionPolytope{1})
  true
end

function (==)(a::ExtrusionPolytope{0},b::ExtrusionPolytope{0})
  true
end

function get_vertex_coordinates(p::ExtrusionPolytope)
  p.vertex_coords
end

function get_edge_tangent(p::ExtrusionPolytope)
  _edge_tangents(Float64,p.dface)
end

function get_facet_normal(p::ExtrusionPolytope)
  p.face_normals
end

function get_facet_orientations(p::ExtrusionPolytope)
  p.face_orientations
end

function get_vertex_permutations(p::ExtrusionPolytope)
  p.vertex_perms
end

function get_face_vertex_permutations(p::ExtrusionPolytope,d::Integer)
  r = get_dimrange(p,d)
  perms = get_face_vertex_permutations(p)
  perms[r]
end

function get_face_vertex_permutations(p::ExtrusionPolytope)
  p.face_vertex_perms
end

function is_simplex(p::ExtrusionPolytope)
  all(Tuple(p.extrusion) .== TET_AXIS)
end

function is_n_cube(p::ExtrusionPolytope)
  all(Tuple(p.extrusion) .== HEX_AXIS)
end

function is_simplex(p::ExtrusionPolytope{0})
  true
end

function is_n_cube(p::ExtrusionPolytope{0})
  true
end

function is_simplex(p::ExtrusionPolytope{1})
  true
end

function is_n_cube(p::ExtrusionPolytope{1})
  true
end

function simplexify(p::ExtrusionPolytope)
  @notimplemented
end

function simplexify(p::ExtrusionPolytope{0})
  [[1,],], VERTEX
end

function simplexify(p::ExtrusionPolytope{1})
  [[1,2],], SEGMENT
end

function simplexify(p::ExtrusionPolytope{2})
  if p == QUAD
    ([[1,2,3],[2,3,4]], TRI)
  elseif p == TRI
    ([[1,2,3],], TRI)
  else
     @notimplemented
  end
end

function simplexify(p::ExtrusionPolytope{3})
  if p == HEX
    simplices = [
      [1,2,3,7], [1,2,5,7], [2,3,4,7],
      [2,4,7,8], [2,5,6,7], [2,6,7,8]]
    (simplices, TET)
  elseif p == TET
    simplices = [[1,2,3,4],]
    (simplices, TET)
  else
     @notimplemented
  end
end

function Base.show(io::IO,p::ExtrusionPolytope)
  if p == VERTEX
    s = "VERTEX"
  elseif p==SEGMENT
    s = "SEGMENT"
  elseif p==QUAD
    s = "QUAD"
  elseif p==TRI
    s = "TRI"
  elseif p==TET
    s = "TET"
  elseif p== HEX
    s = "HEX"
  elseif p==WEDGE
    s = "WEDGE"
  elseif p==PYRAMID
    s = "PYRAMID"
  else
    s = nameof(typeof(p))
  end
  print(io,s)
end

"""
    get_extrusion(p::ExtrusionPolytope)

Equivalent to `p.extrusion`.
"""
get_extrusion(p::ExtrusionPolytope) = p.extrusion



# Helpers for the construction of a DFace

function DFace(extrusion::Point{D,Int}) where D
  zerop = zero(Point{D,Int})
  pol_nfs, pol_dim = _polytopenfaces(zerop, extrusion)
  nfs_id = Dict{NFace{D},Int}(nf => i for (i, nf) in enumerate(pol_nfs))
  nf_nfs, nf_dimranges = _polytopemesh(pol_nfs, nfs_id)
  dims = _nface_to_nfacedim(length(pol_nfs), pol_dim)
  nf_dims = _nf_dims(nf_nfs,nf_dimranges)
  DFace{D}(extrusion, pol_nfs, pol_dim, dims, nf_nfs, nf_dimranges, nf_dims)
end

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

  D = num_components(extrusion)
  zerop = zero(Point{D,Int})
  nf_nfs = Vector{NFace{D}}(undef,0)
  _nfaceboundary!(anchor, zerop, extrusion, true, nf_nfs)

  for i in 1:D
    sort!(nf_nfs, by = x -> x.anchor[i])
  end
  for i in 1:D
    sort!(nf_nfs, by = x -> x.extrusion[i])
  end
  sort!(nf_nfs, by = x -> count(p->p!=0,x.extrusion))

  numnfs = length(nf_nfs)
  nfsdim = [_nfdim(nf_nfs[i].extrusion) for i = 1:numnfs]
  dnf = _nfdim(extrusion)
  dimnfs = Array{UnitRange{Int},1}(undef, dnf + 1)
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

  D = num_components(extrusion)
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
        list = _nfaceboundary!(tetp + edim, extrusion, newext, false, list)
      end
      list = _nfaceboundary!( anchor, extrusion + edim * curex, newext, false, list)

    end
  end

  list
end

function _newext(newext,i)
  m = zero(Mutable(newext))
  D = num_components(newext)
  for j in 1:D
    m[j] = j == i ? 0 : newext[j]
  end
  Point(m)
end

function _edim(newext,i)
  m = zero(Mutable(newext))
  D = num_components(newext)
  for j in 1:D
    m[j] = j == i ? 1 : 0
  end
  Point(m)
end

function _tetp(anchor,i)
  m = zero(Mutable(anchor))
  D = num_components(anchor)
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

function _precompute_m_n_to_mface_to_nface(p::DFace{D}) where D
  m_n_to_mface_to_nface = Matrix{Vector{Vector{Int}}}(undef,D+1,D+1)

  for m in 0:D
    for n in 0:D
      if n <= m
       mface_to_nface = _dimfrom_fs_dimto_fs(p,m,n)
      else
       mface_to_nface = _dimfrom_fs_dimto_fs_dual(p,m,n)
      end
      m_n_to_mface_to_nface[m+1,n+1] = mface_to_nface
    end
  end

  m_n_to_mface_to_nface

end

function _dimfrom_fs_dimto_fs_dual(p::DFace,dimfrom,dimto)
  @assert dimfrom < dimto
  tface_to_ffaces = _dimfrom_fs_dimto_fs(p,dimto,dimfrom)
  nffaces = _num_nfaces(p,dimfrom)
  fface_to_tfaces = [Int[] for in in 1:nffaces]
  for (tface,ffaces) in enumerate(tface_to_ffaces)
    for fface in ffaces
      push!(fface_to_tfaces[fface],tface)
    end
  end
  fface_to_tfaces
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
  b = zero(Mutable(Point{d,Int}))
  D = num_components(a)
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
    vc = Tuple(vx.anchor)
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


function _nfaces_vertices(::Type{T},p::DFace,d::Integer) where T
  nc = _num_nfaces(p,d)
  verts = _vertices_coordinates(T,p)
  faces_vs = _dimfrom_fs_dimto_fs(p,d,0)
  cfvs = collect(lazy_map(Broadcasting(Reindex(verts)),faces_vs))
end

# Return the n-faces vertices coordinates array for a given n-face dimension
function _nfaces_vertices(::Type{T},p::ExtrusionPolytope,d::Integer) where T
  _nfaces_vertices(T,p.dface,d)
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
      for d in 1:num_components(vi)
        v[i-1,d] = vi[d]
      end
    end
    _n = nullspace(v)
    n = Point{D,T}(_n)
    n = n * 1 / sqrt(dot(n, n))
    ext_v = _vertex_not_in_facet(p, i_f, nf_vs)
    v3 = vs[nf_vs[i_f][1]] - vs[ext_v]
    f_or = 1
    if dot(v3, n) < 0.0
      n *= -1
      f_or = -1
    end
  elseif (length(p.extrusion) == 1)
    ext_v = _vertex_not_in_facet(p, i_f, nf_vs)
    n = vs[nf_vs[i_f][1]] - vs[ext_v]
    n = n .* 1 / dot(n, n)
    f_or = 1
  else
    @unreachable "O-dim polytopes do not have properly define outward facet normals"
  end
  n, f_or
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

function _admissible_face_vertex_permutations(p::DFace{D}) where D
  perms = Vector{Vector{Int}}[]
  _admissible_face_vertex_permutations_fill!(perms,p,Val{0}())
  perms
end

function  _admissible_face_vertex_permutations_fill!(perms,p::DFace{D},::Val{d}) where {D,d}
  faceids = p.dimranges[d+1]
  for faceid in faceids
    f = DFace{d}(p,faceid)
    f_perms = _precompute_vertex_perms_if_possible(f)
    push!(perms,f_perms)
  end
  _admissible_face_vertex_permutations_fill!(perms,p,Val{d+1}())
  nothing
end

function  _admissible_face_vertex_permutations_fill!(perms,p::DFace{D},::Val{D}) where D
  f_perms = _precompute_vertex_perms_if_possible(p)
  push!(perms,f_perms)
  nothing
end

function _precompute_vertex_perms_if_possible(p::DFace{D}) where D
  if D in (0,1)
    perms = _admissible_permutations_simplex(p)
  elseif D == 2
    perms = _admissible_permutations(p)
  else
    perms = _admissible_permutations_identity(p)
  end
  perms
end

# It generates all the admissible permutations of nodes that lead to an
# admissible polytope
function _admissible_permutations(p::DFace{D}) where D
  if D > 3
    @warn "Computing permutations for a polytope of dim > 3 is overkill"
  end
  if D in (0,1) || all( map(i->i==TET_AXIS,Tuple(p.extrusion)[2:end]) )
    perms = _admissible_permutations_simplex(p)
  elseif all( map(i->i==HEX_AXIS,Tuple(p.extrusion)[2:end]))
    perms = _admissible_permutations_n_cube(p)
  else
    @notimplemented "admissible vertex permutations only implemented for simplices and n-cubes"
  end
  perms
end

function _admissible_permutations_identity(p::DFace{D}) where D
  vs = first(_dimfrom_fs_dimto_fs(p, D, 0))
  num_vs = length(vs)
  l = [i for i = 1:num_vs]
  perms = [l,]
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
  vol = -1
  for perm in perms
    for (j,cj) in enumerate(perm)
      permuted_vertices[j] = vertices[cj]
    end
    jac = _setup_aux_jacobian(grads,permuted_vertices)
    vol_i = convert(Int,abs(det(jac)))
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
  m = zero(Mutable(Point{D,Int}))
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
  m = zero(Mutable(outer(p0,p0)))
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

# Some particular cases

"""
    const VERTEX = Polytope()
"""
const VERTEX = Polytope()

"""
    const SEGMENT = Polytope(HEX_AXIS)
"""
const SEGMENT = Polytope(HEX_AXIS)

# TODO use larger names
"""
    const TRI = Polytope(TET_AXIS,TET_AXIS)
"""
const TRI = Polytope(TET_AXIS,TET_AXIS)

"""
    const QUAD = Polytope(HEX_AXIS,HEX_AXIS)
"""
const QUAD = Polytope(HEX_AXIS,HEX_AXIS)

"""
    const TET = Polytope(TET_AXIS,TET_AXIS,TET_AXIS)
"""
const TET = Polytope(TET_AXIS,TET_AXIS,TET_AXIS)

"""
    const HEX = Polytope(HEX_AXIS,HEX_AXIS,HEX_AXIS)
"""
const HEX = Polytope(HEX_AXIS,HEX_AXIS,HEX_AXIS)

"""
    const WEDGE = Polytope(TET_AXIS,TET_AXIS,HEX_AXIS)
"""
const WEDGE = Polytope(TET_AXIS,TET_AXIS,HEX_AXIS)

"""
    const PYRAMID = Polytope(HEX_AXIS,HEX_AXIS,TET_AXIS)
"""
const PYRAMID = Polytope(HEX_AXIS,HEX_AXIS,TET_AXIS)
