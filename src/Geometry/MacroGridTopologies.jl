
"""
    struct MacroGridTopology{Dc,Dp} <: GridTopology{Dc,Dp}

Topology of the coarse (macro) faces of a partition of the cells of a fine `GridTopology`.
Every fine d-face belongs to exactly one coarse face, of dimension `>= d`.

The function `get_faces(mt,n,m)` with `n > m` gives the coarse m-faces in the closure of 
the fine n-faces of each coarse n-face. 

## Constructors:

    MacroGridTopology(fine_topo, cell_to_ccell = 1:num_cells(fine_topo); facet_to_tag = zeros)

where `facet_to_tag` is a label (`>= 0`) per fine facet. It is used on the boundary to split 
it into several coarse facets. With the default, connected boundary facets of a coarse cell 
form a single coarse facet, so domain corners (in 3D, also edges) are not coarse vertices (edges).

"""
struct MacroGridTopology{Dc,Dp,A<:GridTopology{Dc,Dp}} <: GridTopology{Dc,Dp}
  fine_topo :: A
  d_to_dface_to_cdim  :: Vector{Vector{Int8}}
  d_to_dface_to_cface :: Vector{Vector{Int32}}
  d_to_ncfaces :: Vector{Int}
  n_m_to_nface_to_mfaces :: Matrix{Table{Int32,Vector{Int32},Vector{Int32}}}
  n_d_to_cface_to_dfaces :: Matrix{Table{Int32,Vector{Int32},Vector{Int32}}}
end

function MacroGridTopology(
  fine_topo::GridTopology{Dc,Dp},
  cell_to_ccell::AbstractVector{<:Integer} = Base.OneTo(Int32(num_cells(fine_topo)));
  facet_to_tag::AbstractVector{<:Integer} = Fill(Int32(0),num_faces(fine_topo,Dc-1))
) where {Dc,Dp}
  facet_to_cfacet = _generate_coarse_facets(fine_topo,cell_to_ccell,facet_to_tag)
  d_to_dface_to_cdim, d_to_dface_to_cface = _generate_coarse_faces(fine_topo,cell_to_ccell,facet_to_cfacet)
  # Every coarse d-face contains a fine d-face
  d_to_ncfaces = [
    Int(maximum(_cface_or_zero.(d_to_dface_to_cdim[d+1],d_to_dface_to_cface[d+1],d);init=Int32(0)))
    for d in 0:Dc
  ]
  T = Table{Int32,Vector{Int32},Vector{Int32}}
  return MacroGridTopology{Dc,Dp,typeof(fine_topo)}(
    fine_topo, d_to_dface_to_cdim, d_to_dface_to_cface, d_to_ncfaces,
    Matrix{T}(undef,Dc+1,Dc+1), Matrix{T}(undef,Dc+1,Dc+1)
  )
end

# GridTopology API

num_faces(mt::MacroGridTopology,d::Integer) = mt.d_to_ncfaces[d+1]

function get_faces(mt::MacroGridTopology,n::Integer,m::Integer)
  _setup_coarse_faces!(mt,n,m)
  mt.n_m_to_nface_to_mfaces[n+1,m+1]
end

function get_vertex_coordinates(mt::MacroGridTopology)
  coords = get_vertex_coordinates(mt.fine_topo)
  lazy_map(Reindex(coords),get_coarse_vertex_to_vertex(mt))
end

function get_polytopes(::MacroGridTopology)
  @notimplemented "The faces of a MacroGridTopology are not polytopes"
end

function get_cell_type(::MacroGridTopology)
  @notimplemented "The faces of a MacroGridTopology are not polytopes"
end

function compute_isboundary_face(mt::MacroGridTopology{Dc},d::Integer) where Dc
  d == Dc && return fill(false,num_faces(mt,Dc))
  facet_to_isboundary = generate_facet_to_isboundary(get_faces(mt,Dc-1,Dc))
  d == Dc-1 && return facet_to_isboundary
  generate_face_to_isboundary(facet_to_isboundary,get_faces(mt,d,Dc-1))
end

# Fine <-> coarse

get_grid_topology(mt::MacroGridTopology) = mt.fine_topo

"""
    get_fine_to_coarse_faces(mt::MacroGridTopology,d) -> (dface_to_cdim, dface_to_cface)

Dimension and id of the coarse face of each fine `d`-face.
"""
function get_fine_to_coarse_faces(mt::MacroGridTopology,d::Integer)
  (mt.d_to_dface_to_cdim[d+1], mt.d_to_dface_to_cface[d+1])
end

"""
    get_coarse_to_fine_faces(mt::MacroGridTopology,n,d) -> cface_to_dfaces

Fine `d`-faces (`d <= n`) of each coarse `n`-face.
"""
function get_coarse_to_fine_faces(mt::MacroGridTopology,n::Integer,d::Integer)
  @check d <= n
  if !isassigned(mt.n_d_to_cface_to_dfaces,n+1,d+1)
    mt.n_d_to_cface_to_dfaces[n+1,d+1] = _coarse_to_fine_faces(
      mt.d_to_dface_to_cdim[d+1],mt.d_to_dface_to_cface[d+1],n,num_faces(mt,n)
    )
  end
  return mt.n_d_to_cface_to_dfaces[n+1,d+1]
end

"""
    get_coarse_vertex_to_vertex(mt::MacroGridTopology) -> cvertex_to_vertex

Fine vertex of each coarse vertex.
"""
function get_coarse_vertex_to_vertex(mt::MacroGridTopology)
  cvertex_to_vertices = get_coarse_to_fine_faces(mt,0,0)
  @check length(cvertex_to_vertices.data) == length(cvertex_to_vertices)
  return cvertex_to_vertices.data
end

# Construction

function _setup_coarse_faces!(mt::MacroGridTopology,n,m)
  isassigned(mt.n_m_to_nface_to_mfaces,n+1,m+1) && return nothing
  if n == m
    mt.n_m_to_nface_to_mfaces[n+1,n+1] = identity_table(Int32,Int32,num_faces(mt,n))
  elseif n < m
    _setup_coarse_faces!(mt,m,n)
  else
    nface_to_mfaces = _generate_coarse_incidence(
      mt.fine_topo,n,m,get_coarse_to_fine_faces(mt,n,n),get_fine_to_coarse_faces(mt,m)...
    )
    mt.n_m_to_nface_to_mfaces[n+1,m+1] = nface_to_mfaces
    mt.n_m_to_nface_to_mfaces[m+1,n+1] = Arrays.inverse_table(nface_to_mfaces,num_faces(mt,m))
  end
  return nothing
end

# Fine facet -> coarse facet (0 if interior to a coarse cell)
function _generate_coarse_facets(topo::GridTopology,cell_to_ccell,facet_to_tag)
  D = num_cell_dims(topo)
  facet_to_cells = get_faces(topo,D-1,D)
  nfacets = length(facet_to_cells)

  # Key: (c1,c2), c1 < c2, between two coarse cells; (c,-tag) on the boundary; (0,0) inside
  facet_to_key = Vector{Tuple{Int32,Int32}}(undef,nfacets)
  for f in 1:nfacets
    cells = dataview(facet_to_cells,f)
    c1 = Int32(cell_to_ccell[cells[1]])
    facet_to_key[f] = if length(cells) == 1
      (c1, -Int32(facet_to_tag[f]))
    else
      @check length(cells) == 2
      c2 = Int32(cell_to_ccell[cells[2]])
      c1 == c2 ? (Int32(0),Int32(0)) : minmax(c1,c2)
    end
  end
  isactive(f) = !iszero(first(facet_to_key[f]))
  samekey(f,g) = isactive(f) && facet_to_key[f] == facet_to_key[g]

  sets = IntDisjointSets{Int32}(nfacets)
  if D >= 2
    _union_neighbours!(sets,0,get_faces(topo,D-2,D-1),samekey)
  end
  facet_to_cfacet, _ = _number_classes(sets,isactive,1:nfacets)
  return facet_to_cfacet
end

# For every fine d-face, d < D: sorted coarse facets containing a fine facet incident to it
function _generate_incident_coarse_facets(topo::GridTopology,facet_to_cfacet)
  D = num_cell_dims(topo)
  ncfacets = Int(maximum(facet_to_cfacet;init=Int32(0)))
  facet_to_K = Table(facet_to_cfacet,collect(Int32(1):Int32(length(facet_to_cfacet)+1)))
  return [Arrays.merge_entries(facet_to_K,get_faces(topo,d,D-1),ncfacets) for d in 0:D-1]
end

# Fine d-face -> (dim, id) of its coarse face, d = 0:D
function _generate_coarse_faces(topo::GridTopology,cell_to_ccell,facet_to_cfacet)
  D = num_cell_dims(topo)
  d_to_dface_to_K = _generate_incident_coarse_facets(topo,facet_to_cfacet)

  # Union-find on the fine faces of dimension < D, with global ids offsets[d+1] + dface
  offsets = cumsum([0; [num_faces(topo,d) for d in 0:D-1]])
  sets = IntDisjointSets{Int32}(offsets[end])
  isactive = falses(offsets[end]) # belongs to a coarse face of dimension < D
  for d in 0:D-1
    K, o = d_to_dface_to_K[d+1], offsets[d+1]
    hasK(f) = !isempty(dataview(K,f))
    sameK(f,g) = hasK(f) && dataview(K,f) == dataview(K,g)
    dface_to_cells = get_faces(topo,d,D)
    for f in eachindex(K)
      isactive[o+f] = hasK(f) || _touches_several(dataview(dface_to_cells,f),cell_to_ccell)
    end
    if d < D-1 # incident (d+1)-faces with the same K
      Kup, oup = d_to_dface_to_K[d+2], offsets[d+2]
      dface_to_upfaces = get_faces(topo,d,d+1)
      for f in eachindex(K), g in dataview(dface_to_upfaces,f)
        (hasK(f) && dataview(K,f) == dataview(Kup,g)) && union!(sets,Int32(o+f),Int32(oup+g))
      end
    end
    if d >= 1 # d-faces with the same K sharing a (d-1)-face
      _union_neighbours!(sets,o,get_faces(topo,d-1,d),sameK)
    end
  end

  # Coarse dimension of each class, then dimension-local numbering
  gface_to_dim = vcat([fill(Int8(d),num_faces(topo,d)) for d in 0:D-1]...)
  gface_to_cface, gface_to_cdim = _number_classes(
    sets,i->isactive[i],1:offsets[end];class_dim=gface_to_dim
  )

  d_to_dface_to_cdim  = [fill(Int8(D),num_faces(topo,d)) for d in 0:D]
  d_to_dface_to_cface = [zeros(Int32,num_faces(topo,d)) for d in 0:D]
  for d in 0:D-1
    dface_to_cells = get_faces(topo,d,D)
    for f in 1:num_faces(topo,d)
      i = offsets[d+1] + f
      if isactive[i]
        d_to_dface_to_cdim[d+1][f]  = gface_to_cdim[i]
        d_to_dface_to_cface[d+1][f] = gface_to_cface[i]
      else # interior to a coarse cell
        d_to_dface_to_cface[d+1][f] = cell_to_ccell[first(dataview(dface_to_cells,f))]
      end
    end
  end
  d_to_dface_to_cface[D+1] .= cell_to_ccell
  return d_to_dface_to_cdim, d_to_dface_to_cface
end

# Coarse n-face -> coarse m-faces in the closure of its fine n-faces (n > m)
function _generate_coarse_incidence(topo::GridTopology,n,m,cface_to_nfaces,mface_to_cdim,mface_to_cface)
  @check n > m
  mface_to_cmface = _cface_or_zero.(mface_to_cdim,mface_to_cface,m)
  ncmfaces = Int(maximum(mface_to_cmface;init=Int32(0)))
  # 0 for fine m-faces that are not coarse m-faces, dropped by the merge
  nface_to_mfaces = get_faces(topo,n,m)
  nface_to_cmfaces = Table(mface_to_cmface[nface_to_mfaces.data],nface_to_mfaces.ptrs)
  return Arrays.merge_entries(nface_to_cmfaces,cface_to_nfaces,ncmfaces)
end

_cface_or_zero(cdim,cface,n) = cdim == n ? cface : zero(cface)

function _coarse_to_fine_faces(dface_to_cdim,dface_to_cface,n,ncfaces)
  dface_to_cnface = _cface_or_zero.(dface_to_cdim,dface_to_cface,n)
  data, ptrs = Arrays.inverse_table(dface_to_cnface,Int32(1):Int32(length(dface_to_cnface)+1),ncfaces)
  return Table(data,ptrs)
end

# Union the d-faces f, g (elements offset+f, offset+g) sharing a (d-1)-face with issame(f,g)
function _union_neighbours!(sets,offset,lface_to_dfaces::Table,issame)
  for l in eachindex(lface_to_dfaces)
    dfaces = dataview(lface_to_dfaces,l)
    for j in 2:length(dfaces), i in 1:j-1
      fi, fj = dfaces[i], dfaces[j]
      issame(fi,fj) && union!(sets,Int32(offset+fi),Int32(offset+fj))
    end
  end
end

# Number the classes of the active elements in order of first appearance (0 if inactive).
# With class_dim, numbering is per dimension; the dimension of a class is the largest of its elements.
function _number_classes(sets,isactive,ids;class_dim=nothing)
  n = length(sets)
  root_to_dim = isnothing(class_dim) ? nothing : fill(Int8(-1),n)
  if !isnothing(class_dim)
    for i in ids
      isactive(i) || continue
      r = find_root!(sets,Int32(i))
      root_to_dim[r] = max(root_to_dim[r],class_dim[i])
    end
  end
  root_to_class = zeros(Int32,n)
  dim_to_nclasses = zeros(Int32,isnothing(class_dim) ? 1 : maximum(class_dim;init=Int8(0))+1)
  id_to_class = zeros(Int32,n)
  id_to_dim = fill(Int8(-1),n)
  for i in ids
    isactive(i) || continue
    r = find_root!(sets,Int32(i))
    k = isnothing(class_dim) ? 1 : root_to_dim[r] + 1
    if iszero(root_to_class[r])
      dim_to_nclasses[k] += 1
      root_to_class[r] = dim_to_nclasses[k]
    end
    id_to_class[i] = root_to_class[r]
    id_to_dim[i] = k - 1
  end
  return id_to_class, id_to_dim
end

function _touches_several(cells,cell_to_ccell)
  c = cell_to_ccell[first(cells)]
  return any(x -> cell_to_ccell[x] != c, cells)
end
