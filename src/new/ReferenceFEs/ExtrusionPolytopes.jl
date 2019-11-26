
# Helper type
# n-face of the polytope, i.e., any polytope of lower dimension `N` representing
# its boundary and the polytope itself (for `N` equal to the space dimension `D`)
struct NFace{D}
  anchor::Point{D,Int}
  extrusion::Point{D,Int}
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
  nfaces::Vector{NFace}
  nfacesdim::Vector{UnitRange{Int}}
  nf_nfs::Vector{Vector{Int}}
  nf_dim::Vector{Vector{UnitRange{Int}}}
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

function ExtrusionPolytope(extrusion::NTuple{N,Int}) where N
  ExtrusionPolytope(Point{N,Int}(extrusion))
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

"""
Constant to be used in order to indicate a hex-like extrusion axis.
"""
const HEX_AXIS = 1

"""
Constant to be used in order to indicate a tet-like extrusion axis.
"""
const TET_AXIS = 2

# Implementation of the interface

get_faces(p::ExtrusionPolytope) = p.nf_nfs

get_dimranges(p::ExtrusionPolytope) = p.nfacesdim

function Polytope{0}(p::ExtrusionPolytope,Dfaceid::Integer)
  @assert Dfaceid in 1:num_vertices(p)
  _vertex()
end

function Polytope{D}(p::ExtrusionPolytope,Dfaceid::Integer)::ExtrusionPolytope{D} where D
  reffaces = _nface_ref_polytopes(p)
  offset = get_offset(p,D)
  id = Dfaceid + offset
  reffaces[id]
end

function Polytope{D}(p::ExtrusionPolytope{D},Dfaceid::Integer) where D
  @assert Dfaceid == 1
  p
end

function (==)(a::ExtrusionPolytope{D},b::ExtrusionPolytope{D}) where D
  a.extrusion == b.extrusion
end

function get_vertex_coordinates(p::ExtrusionPolytope)
  _vertices_coordinates(p)
end

function get_edge_tangents(p::ExtrusionPolytope)
  _edge_tangents(p)
end

function get_facet_normals(p::ExtrusionPolytope)
  n, _ = _face_normals(p)
  n
end

function get_facet_orientations(p::ExtrusionPolytope)
  _, or = _face_normals(p)
  or
end

function get_vertex_permutations(p::ExtrusionPolytope)
  _admissible_permutations(p)
end

function is_simplex(p::ExtrusionPolytope)
  all(p.extrusion.array .== TET_AXIS)
end

function is_n_cube(p::ExtrusionPolytope)
  all(p.extrusion.array .== HEX_AXIS)
end

function is_simplex(p::ExtrusionPolytope{0})
  true
end

function is_n_cube(p::ExtrusionPolytope{0})
  true
end

# Helpers

function ExtrusionPolytope(extrusion::Point{D,Int}) where D
  zerop = Point{D,Int}(zeros(Int64, D))
  pol_nfs_dim = _polytopenfaces(zerop, extrusion)
  pol_nfs = pol_nfs_dim[1]
  pol_dim = pol_nfs_dim[2]
  nfs_id = Dict(nf => i for (i, nf) in enumerate(pol_nfs))
  nf_nfs_dim = _polytopemesh(pol_nfs, nfs_id)
  nf_nfs = nf_nfs_dim[1]
  nf_dim = nf_nfs_dim[2]
  ExtrusionPolytope{D}(extrusion, pol_nfs, pol_dim, nf_nfs, nf_dim)
end

function ExtrusionPolytope(extrusion::Point{0,Int})
  _vertex()
end

# Generates the array of n-faces of a polytope
function _polytopenfaces(anchor::Point{D,Int}, extrusion::Point{D,Int}) where D
  dnf = _nfdim(extrusion)
  zerop = Point{D,Int}(zeros(Int64, D))
  nf_nfs = []
  nf_nfs = _nfaceboundary!(anchor, zerop, extrusion, true, nf_nfs)
  [sort!(nf_nfs, by = x -> x.anchor[i]) for i = 1:length(extrusion)]
  [sort!(nf_nfs, by = x -> x.extrusion[i]) for i = 1:length(extrusion)]
  [sort!(nf_nfs, by = x -> sum(x.extrusion))]
  numnfs = length(nf_nfs)
  nfsdim = [_nfdim(nf_nfs[i].extrusion) for i = 1:numnfs]
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
  return [nf_nfs, dimnfs]
end

_nfdim(a::Point{D,Int}) where D = sum([a[i] > 0 ? 1 : 0 for i = 1:D])

# Generates the list of n-face of a polytope the d-faces for 0 <= d <n on its
# boundary
function _nfaceboundary!(
  anchor::Point{D,Int},
  extrusion::Point{D,Int},
  extend::Point{D,Int},
  isanchor::Bool,
  list) where D

  newext = extend
  list = [list..., NFace{D}(anchor, extrusion)]
  for i = 1:D
    curex = newext[i]
    if (curex > 0) # Perform extension
      func1 = (j -> j == i ? 0 : newext[j])
      newext = Point{D,Int}([func1(i) for i = 1:D])
      func2 = (j -> j == i ? 1 : 0)
      edim = Point{D,Int}([func2(i) for i = 1:D])
      func3 = (j -> j >= i ? anchor[j] : 0)
      tetp = Point{D,Int}([func3(i) for i = 1:D]) + edim
      if (curex == 1) # Quad extension
        list = _nfaceboundary!(anchor + edim, extrusion, newext, false, list)
      elseif (isanchor)
        list = _nfaceboundary!(tetp, extrusion, newext, false, list)
      end
      list = _nfaceboundary!(
        anchor,
        extrusion + edim * curex,
        newext,
        false,
        list)
    end
  end
  return list
end

# Provides for all n-faces of a polytope the d-faces for 0 <= d <n on its
# boundary (e.g., given a face, it provides the set of edges and corners on its
# boundary) using the global n-face numbering of the base polytope
function _polytopemesh(nfaces::Vector{NFace{D}}, nfaceid::Dict) where D
  num_nfs = length(nfaces)
  nfnfs = Vector{Vector{Int64}}(undef, num_nfs)
  nfnfs_dim = Vector{Vector{UnitRange{Int64}}}(undef, num_nfs)
  for (inf, nf) in enumerate(nfaces)
    nfs_dim_nf = _polytopenfaces(nf.anchor, nf.extrusion)
    nf_nfs = nfs_dim_nf[1]
    dimnfs = nfs_dim_nf[2]
    nfnfs[inf] = [get(nfaceid, nf, nf) for nf in nf_nfs]
    nfnfs_dim[inf] = dimnfs
  end
  return [nfnfs, nfnfs_dim]
end

# Returns an array with the reference polytopes for all n-faces (undef for vertices)
function _nface_ref_polytopes(p::ExtrusionPolytope)
  function _eliminate_zeros(a)
    b = Int[]
    for m in a
      if (m != 0)
        push!(b, m)
      end
    end
    return Tuple(b)
  end
  nf_ref_p = Vector{ExtrusionPolytope}(undef, length(p.nfaces))
  ref_nf_ps = ExtrusionPolytope[]
  v = _vertex()
  for (i_nf, nf) in enumerate(p.nfaces)
    r_ext = _eliminate_zeros(nf.extrusion)
    if r_ext != ()
      k = 0
      for (i_p, ref_p) in enumerate(ref_nf_ps)
        if r_ext == ref_p.extrusion
          k = i_p
          nf_ref_p[i_nf] = ref_p
        end
      end
      if k == 0
        ref_p = ExtrusionPolytope(r_ext)
        push!(ref_nf_ps, ref_p)
        k = length(ref_nf_ps) + 1
        nf_ref_p[i_nf] = ref_p
      end
    else
        nf_ref_p[i_nf] = v
    end
  end
  return nf_ref_p
end

function _vertex()
  ext = ()
  nfdim = [[1:1]]
  nfnfs = [[1]]
  nfanc = Point{0,Int}()
  nf = NFace{0}(nfanc,nfanc)
  nfs = [nf]
  return ExtrusionPolytope{0}(ext, nfs, [1:1], nfnfs, nfdim)
end

function _dimfrom_fs_dimto_fs(p::ExtrusionPolytope, dim_from::Int, dim_to::Int)
  @assert dim_to <= dim_from
  dim_from += 1
  dim_to += 1
  dffs_r = p.nf_dim[end][dim_from]
  dffs_dtfs = Vector{Vector{Int}}(undef, dffs_r[end] - dffs_r[1] + 1)
  offs = p.nf_dim[end][dim_to][1] - 1
  for (i_dff, dff) in enumerate(dffs_r)
    dff_nfs = p.nf_nfs[dff]
    dff_dtfs_r = p.nf_dim[dff][dim_to]
    dff_dtfs = dff_nfs[dff_dtfs_r]
    dffs_dtfs[i_dff] = dff_dtfs .- offs
    # @santiagobadia : With or without offset ?
  end
  return dffs_dtfs
end

# It generates the list of coordinates of all vertices in the polytope. It is
# assumed that the polytope has the bounding box [0,1]**dim
function _vertices_coordinates(p::ExtrusionPolytope{D}) where D
  vs = _dimfrom_fs_dimto_fs(p, D, 0)[1]
  vcs = Point{D,Float64}[]
  for i = 1:length(vs)
    cs = convert(Vector{Float64}, [p.nfaces[vs[i]].anchor...])
    push!(vcs, cs)
  end
  return vcs
end

# Return the n-faces vertices coordinates array for a given n-face dimension
function _nfaces_vertices(p,d)
  nc = _num_nfaces(p,d)
  verts = _vertices_coordinates(p)
  faces_vs = _dimfrom_fs_dimto_fs(p,d,0)
  cfvs = collect(LocalToGlobalArray(faces_vs,verts))
end

# Returns number of nfaces of dimension dim
function _num_nfaces(polytope::ExtrusionPolytope, dim::Integer)
  n = length(polytope.nf_dim)
  k = 0
  for nface = 1:n
    d = length(polytope.nf_dim[nface]) - 1
    if d == dim
      k += 1
    end
  end
  k
end

# It generates the outwards normals of the facets of a polytope. It returns two
# arrays, the first one being the outward normal and the second one the orientation.
function _face_normals(p::ExtrusionPolytope{D}) where D
  nf_vs = _dimfrom_fs_dimto_fs(p, D - 1, 0)
  vs = _vertices_coordinates(p)
  f_ns = Point{D,Float64}[]
  f_os = Int[]
  for i_f = 1:length(p.nf_dim[end][end-1])
    n, f_o = _facet_normal(p, nf_vs, vs, i_f)
    push!(f_ns, Point{D,Float64}(n))
    push!(f_os, f_o)
  end
  return f_ns, f_os
end

function _face_normals(p::ExtrusionPolytope{0})
  (VectorValue{0,Float64}[], Int[])
end

function _face_normals(p::ExtrusionPolytope{1})
  (VectorValue{1,Float64}[(-1),(1)], [1,1])
end

# It generates the tangent vectors for polytope edges.
function _edge_tangents(p::ExtrusionPolytope{D}) where D
  ed_vs = _nfaces_vertices(p,1)
  return ts = [(t = vs[2]-vs[1])/norm(t) for vs in ed_vs ]
end

function _edge_tangents(p::ExtrusionPolytope{0})
  VectorValue{0,Float64}[]
end

function _facet_normal(p::ExtrusionPolytope{D}, nf_vs, vs, i_f) where D
  if (length(p.extrusion) > 1)
    v = Float64[]
    for i = 2:length(nf_vs[i_f])
      vi = vs[nf_vs[i_f][i]] - vs[nf_vs[i_f][1]]
      push!(v, vi...)
    end
    n = nullspace(transpose(reshape(v, D, length(nf_vs[i_f]) - 1)))
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
    n = vs[nf_vs[i_f][1]] - vs[ext_v]
    n = n .* 1 / dot(n, n)
    f_or = 1
  else
    error("O-dim polytopes do not have properly define outward facet normals")
  end
  return n, f_or
end

function _vertex_not_in_facet(p, i_f, nf_vs)
  for i in p.nf_dim[end][1]
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
end


# It generates all the admissible permutations of nodes that lead to an
# admissible polytope
function _admissible_permutations(p::ExtrusionPolytope)
  p_dims = length(p.extrusion)
  p_vs = _dimfrom_fs_dimto_fs(p, p_dims, 0)
  vs = p.nfaces[p_vs...]
  num_vs = length(vs)
  ext = p.extrusion
  l = [i for i = 1:num_vs]
  permuted_polytopes = Vector{Int}[]
  for c in Combinatorics.permutations(l)#, p_dims + 1)
    admissible_polytope = true
    c1 = vs[c[1]].anchor
    for j = 2:p_dims+1
      c2 = vs[c[j]].anchor
      if (!_are_nodes_connected(c1, c2, ext))
        admissible_polytope = false
      end
    end
    if (admissible_polytope)
      push!(permuted_polytopes, c)
    end
  end
  return permuted_polytopes
end

# Auxiliary function that determines whether two nodes are connected
function _are_nodes_connected(c1, c2, ext)
  sp_dims = length(c1)
  d = zeros(length(c1))
  for i = 1:length(d)
    d[i] = c2[i] - c1[i]
  end
  dn = sum(d .* d)
  connected = false
  if (dn == 1)
    connected = true
  else
    k = 0
    for j = 1:sp_dims
      if (ext[j] == 2)
        if (c1[j] == 1 || c2[j] == 1)
          k = j
        end
      end
    end
    for l = 1:k
      d[l] = 0
    end
    dn = sum(d .* d)
    if (dn == 0)
      connected = true
    end
  end
  return connected
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

