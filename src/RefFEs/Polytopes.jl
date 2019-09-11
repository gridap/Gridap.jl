module Polytopes

using Gridap
using Gridap.Helpers
using StaticArrays
using Base.Cartesian
using LinearAlgebra

using Combinatorics

export Polytope
export NFace
export HEX_AXIS, TET_AXIS
export num_nfaces
export facet_normals

export space_dim
export dim
export anchor
export extrusion
export nfaces
export nfaces_dim
export nface_connections
export nf_nfs
export nf_dim

export num_nfaces
export nface_ref_polytopes
export generate_admissible_permutations
export equidistant_interior_nodes_coordinates
export vertices_coordinates
export facet_normals

# Module constants

const HEX_AXIS = 1
const TET_AXIS = 2

# Concrete implementations

"""
n-face of the polytope, i.e., any polytope of lower dimension `N` representing
its boundary and the polytope itself (for `N` equal to the space dimension `D`)
"""
struct NFace{D}
  anchor::Point{D,Int}
  extrusion::Point{D,Int}
end

"""
Aggregation of all n-faces that compose the polytope boundary and the polytope
itself, the classification of n-faces with respect to their dimension and type
"""
struct Polytope{D}
  extrusion::Point{D,Int}
  nfaces::Vector{NFace}
  nf_nfs::Vector{Vector{Int64}}
  nf_dim::Vector{Vector{UnitRange{Int64}}}
end

"""
Constructs a `Polytope` given the type of extrusion, i.e., HEX_AXIS (=1) for "hex" extrusion
and TET_AXIS (=2) for "tet" extrusion
"""
function Polytope(extrusion::Vararg{Int,N}) where N
  return Polytope(Point{N,Int}(extrusion))
end

function Polytope(extrusion::NTuple{N,Int}) where N
  return Polytope(extrusion...)
end

function Polytope(extrusion::Point{D,Int}) where D
  zerop = Point{D,Int}(zeros(Int64, D))
  pol_nfs_dim = _polytopenfaces(zerop, extrusion)
  pol_nfs = pol_nfs_dim[1]
  pol_dim = pol_nfs_dim[2]
  nfs_id = Dict(nf => i for (i, nf) in enumerate(pol_nfs))
  nf_nfs_dim = _polytopemesh(pol_nfs, nfs_id)
  nf_nfs = nf_nfs_dim[1]
  nf_dim = nf_nfs_dim[2]
  Polytope{D}(extrusion, pol_nfs, nf_nfs, nf_dim)
end

"""
Provides the dimension of the environment space in which the n-face is defined
"""
space_dim(::NFace{D}) where {D} = D

anchor(nf::NFace) = nf.anchor

"""
Provides the dimension of the n-face, i.e., n
"""
dim(nf::NFace) = _nfdim(nf.extrusion)

extrusion(nf::NFace) = nf.extrusion

dim(::Polytope{D}) where {D} = D

extrusion(p::Polytope) = p.extrusion

nfaces(p::Polytope) = p.nfaces

nfaces_dim(p::Polytope, d) = p.nf_dim[end][d+1]

"""
# It provides for every df-face in the polytope all its dt-faces
# We use dim-wise numbering, i.e., we start numbering from 1 at every dim
"""
nface_connections(p::Polytope, dfrom, dto) = _dimfrom_fs_dimto_fs(p, dfrom, dto)

# @santiagobadia : I would prefer not to make public the following ones
"""
Provides the label of the n-faces that are a subset (on the boundary) or
equal to a given n-face. We use a global numbering for n-faces in the returned
array. n-faces are sorted by increasing dimension.
"""
nf_nfs(p::Polytope) = p.nf_nfs

"""
Provides for the result of the `nf_nfs` method, i.e., n-faces in an n-face,
the range of n-faces for every dimension lower or equal than the n-face dimension.
"""
nf_dim(p::Polytope) = p.nf_dim

"""
Returns number of nfaces of dimension dim
"""
function num_nfaces(polytope::Polytope, dim::Integer)
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

"""
Returns number of nfaces
"""
num_nfaces(polytope::Polytope) = length(polytope.nfaces)

"""
# Returns an array with the reference polytopes for all n-faces (undef for vertices)
"""
function nface_ref_polytopes(p::Polytope)
  function _eliminate_zeros(a)
    b = Int[]
    for m in a
      if (m != 0)
        push!(b, m)
      end
    end
    return Tuple(b)
  end
  nf_ref_p = Vector{Polytope}(undef, length(p.nfaces))
  ref_nf_ps = Polytope[]
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
        ref_p = Polytope(r_ext)
        push!(ref_nf_ps, ref_p)
        k = length(ref_nf_ps) + 1
        nf_ref_p[i_nf] = ref_p
      end
    end
  end
  return nf_ref_p
end

"""
It generates all the admissible permutations of nodes that lead to an
admissible polytope
"""
function generate_admissible_permutations(p::Polytope)
  p_dims = length(p.extrusion)
  p_vs = Gridap.Polytopes._dimfrom_fs_dimto_fs(p, p_dims, 0)
  vs = p.nfaces[p_vs...]
  num_vs = length(vs)
  ext = p.extrusion
  l = [i for i = 1:num_vs]
  # @santiagobadia : Here we have to decide how we want this info stored
  permuted_polytopes = Vector{Int}[]
  for c in Combinatorics.permutations(l, p_dims + 1)
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

"""
It generates the set of nodes (its coordinates) in the interior of an n-face,
for a given order. The node coordinates are the ones for a equispace case.
"""
function equidistant_interior_nodes_coordinates(p::Polytope{D}, order::Int) where D
  _order = order*ones(Int64,dim(p))
  equidistant_interior_nodes_coordinates(p, _order)
end

function equidistant_interior_nodes_coordinates(p::Polytope{D}, order::Vector{Int}) where D

  if (all(extrusion(p).array .== HEX_AXIS) || all(order .== order[1]))
    ns = _interior_nodes_int_coords(p, order)
    return ns_float = _interior_nodes_int_to_real_coords(ns, order)
  else
    error("One can consider anisotropic orders on n-cubes only")
  end
end

"""
It generates the list of coordinates of all vertices in the polytope. It is
assumed that the polytope has the bounding box [0,1]**dim
"""
function vertices_coordinates(p::Polytope{D}) where D
  vs = _dimfrom_fs_dimto_fs(p, D, 0)[1]
  vcs = Point{D,Float64}[]
  for i = 1:length(vs)
    cs = convert(Vector{Float64}, [p.nfaces[vs[i]].anchor...])
    push!(vcs, cs)
  end
  return vcs
end

"""
It generates the outwards normals of the facets of a polytope. It returns two
arrays, the first one being the outward normal and the second one the orientation.
"""
function facet_normals(p::Polytope{D}) where D
  nf_vs = _dimfrom_fs_dimto_fs(p, D - 1, 0)
  vs = vertices_coordinates(p)
  f_ns = Point{D,Float64}[]
  f_os = Int[]
  for i_f = 1:length(p.nf_dim[end][end-1])
    n, f_o = _facet_normal(p, nf_vs, vs, i_f)
    push!(f_ns, Point{D,Float64}(n))
    push!(f_os, f_o)
  end
  return f_ns, f_os
end

# Helpers

_nfdim(a::Point{D,Int}) where D = sum([a[i] > 0 ? 1 : 0 for i = 1:D])

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

# Generates the list of n-face of a polytope the d-faces for 0 <= d <n on its
# boundary
function _nfaceboundary!(
  anchor::Point{D,Int},
  extrusion::Point{D,Int},
  extend::Point{D,Int},
  isanchor::Bool,
  list
) where D
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
        list
      )
    end
  end
  return list
end

function _dimfrom_fs_dimto_fs(p::Polytope, dim_from::Int, dim_to::Int)
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

# It generates the set of nodes (its coordinates) in the interior of an n-face,
# for a given order. The node coordinates are `Int` and from 0 to `order` per
# direction
function _interior_nodes_int_coords(p::Polytope{D}, order) where D
  ext = p.extrusion
  _ord = [order...]
  verts = Point{D,Int}[]
  coor = zeros(Int, D)
  _generate_nodes!(D, p.extrusion, _ord, coor, verts)
  return verts
end

# Auxiliary private recursive function to implement _interior_nodes_int_coords
function _generate_nodes!(dim, ext, order, coor, verts)
  ncoo = copy(coor)
  nord = copy(order)
  for i = 1:order[dim]-1
    ncoo[dim] = i
    if dim > 1
      if (ext[dim] == TET_AXIS)
        nord .-= 1
      end
      _generate_nodes!(dim - 1, ext, nord, ncoo, verts)
    else
      push!(verts, Point(ncoo...))
    end
  end
end

# Transforms the int coordinates to float coordinates
function _interior_nodes_int_to_real_coords(nodes, order)
  if length(nodes) > 0
    dim = length(nodes[1])
    cs_float = Point{dim,Float64}[]
    cs = zeros(Float64, dim)
    for cs_int in nodes
      for i = 1:dim
        cs[i] = cs_int[i] / order[i]
      end
      push!(cs_float, Point{dim,Float64}(cs))
    end
  else
    cs_float = Point{0,Float64}[]
  end
  return cs_float
end

function _facet_normal(p::Polytope{D}, nf_vs, vs, i_f) where D
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

end # module Polytopes
