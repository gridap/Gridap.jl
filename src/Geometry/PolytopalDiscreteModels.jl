
"""
    struct PolytopalGridTopology{Dc,Dp,Tp} <: GridTopology{Dc,Dp}
      vertex_coordinates::Vector{Point{Dp,Tp}}
      n_m_to_nface_to_mfaces::Matrix{Table{Int32,Vector{Int32},Vector{Int32}}}
      polytopes::Vector{GeneralPolytope{Dc,Dp,Tp,Nothing}}
    end

Grid topology for polytopal grids.

Constructors: 

    PolytopalGridTopology(topo::UnstructuredGridTopology)
    PolytopalGridTopology(vertex_coordinates,cell_vertices,polytopes)

"""
struct PolytopalGridTopology{Dc,Dp,Tp} <: GridTopology{Dc,Dp}
  vertex_coordinates::Vector{Point{Dp,Tp}}
  n_m_to_nface_to_mfaces::Matrix{Table{Int32,Vector{Int32},Vector{Int32}}}
  polytopes::Vector{GeneralPolytope{Dc,Dp,Tp,Nothing}}
end

function PolytopalGridTopology(
  vertex_coordinates::Vector{<:Point{Dp,T}},
  cell_vertices::Table,
  polytopes::Vector{<:GeneralPolytope{Dc}}
) where {Dc,Dp,T}
  @check length(cell_vertices) == length(polytopes)

  n = Dc+1
  n_vertices = length(vertex_coordinates)
  n_m_to_nface_to_mfaces = Matrix{Table{Int32,Vector{Int32},Vector{Int32}}}(undef,n,n)
  n_m_to_nface_to_mfaces[Dc+1,0+1] = cell_vertices
  vertex_cells = generate_cells_around(cell_vertices,n_vertices)
  n_m_to_nface_to_mfaces[0+1,Dc+1] = vertex_cells

  PolytopalGridTopology{Dc,Dp,T}(
    vertex_coordinates, n_m_to_nface_to_mfaces, polytopes
  )
end

function PolytopalGridTopology(
  vertex_coordinates::Vector{<:Point{Dp,T}},
  d_to_dface_vertices::Vector{<:Table},
  polytopes::Vector{<:GeneralPolytope{Dc}}
) where {Dc,Dp,T}

  n = Dc+1
  n_vertices = length(vertex_coordinates)
  n_m_to_nface_to_mfaces = Matrix{Table{Int32,Vector{Int32},Vector{Int32}}}(undef,n,n)
  for d in 0:Dc
    dface_to_vertices = d_to_dface_vertices[d+1]
    n_m_to_nface_to_mfaces[d+1,0+1] = dface_to_vertices
    vertex_to_dfaces = generate_cells_around(dface_to_vertices,n_vertices)
    n_m_to_nface_to_mfaces[0+1,d+1] = vertex_to_dfaces
  end

  PolytopalGridTopology{Dc,Dp,T}(
    vertex_coordinates, n_m_to_nface_to_mfaces, polytopes
  )
end

function PolytopalGridTopology(topo::UnstructuredGridTopology{Dc}) where Dc
  polys, ctypes = get_polytopes(topo), get_cell_type(topo)
  cell_polys = expand_cell_data(polys,ctypes)
  vertex_coords = get_vertex_coordinates(topo)
  polytopes, cell_nodes = generate_polytopes(
    cell_polys,get_faces(topo,Dc,0),vertex_coords
  )
  PolytopalGridTopology(vertex_coords,cell_nodes,polytopes)
end

#= NOTE: Unstructured -> Polytopal topology conversion

The main idea is that there is a permutation of the local vertex/face numbering 
between ExtrusionPolytopes and GeneralPolytopes. 

In 3D, this is not too bad: Vertex numbering is kept the same, so we can deduce 
the face numbering permutation from the vertex numbering permutation.

In 2D, things are quite more involved... Vertices can potentially undergo two separate 
permutations: First, a permutation that makes them a cyclic graph. Second, a possible 
reversal of the cyclic ordering to ensure positive (counter-clockwise) orientation.

The machinery below hopefully deals with all these issues.
=#

function oriented_cell_nodes(
  cell_nodes::Table{T}, node_coordinates::AbstractVector{<:Point{2}}
) where T
  ptrs = cell_nodes.ptrs
  data = zeros(T, length(cell_nodes.data))
  for cell in 1:length(cell_nodes)
    nodes = dataview(cell_nodes,cell)
    perm  = orient_nodes_positive(view(node_coordinates,nodes); keepfirst=true)
    data[ptrs[cell]:ptrs[cell+1]-1] .= nodes[perm]
  end
  return Table(data,ptrs)
end

function oriented_cell_nodes(
  cell_polys::AbstractArray{<:Polytope{2}}, cell_nodes::Table{T}, node_coordinates::AbstractVector{<:Point{2}}
) where T
  ptrs = cell_nodes.ptrs
  data = zeros(T, length(cell_nodes.data))
  for cell in 1:length(cell_nodes)
    poly  = cell_polys[cell]
    nodes = dataview(cell_nodes,cell)
    pos   = _polytope_is_positive(poly, view(node_coordinates,nodes))
    perm  = _polytope_node_reindex(poly, pos)
    data[ptrs[cell]:ptrs[cell+1]-1] .= nodes[perm]
  end
  return Table(data,ptrs)
end

# Face reindexing between standard and general polytopes
function _polytope_node_reindex(p::Polytope{D}, positive::Bool) where D
  (D == 3) && return collect(1:num_vertices(p))
  if p == TRI
    perm = collect(1:num_vertices(p))
  elseif p == QUAD
    perm = [1,2,4,3]
  else
    @notimplemented
  end
  if !positive
    reverse!(perm)
  end
  return perm
end

function _polytope_is_positive(p::Polytope, coords)
  perm = _polytope_node_reindex(p, true)
  return is_positive_node_orientation(view(coords, perm))
end

function _polytope_face_reindex(p::Polytope{D}, d::Integer, positive::Bool) where D
  n2o_node = _polytope_node_reindex(p, positive)
  iszero(d) && return n2o_node

  q = GeneralPolytope{D}(p)
  old_f2n = get_faces(p,d,0)
  new_f2n = get_faces(q,d,0)

  same_face(f,g) = all(n -> n ∈ g, f)
  n2o_face = zeros(Int,num_faces(p,d))
  for (new, new_face) in enumerate(new_f2n)
    old = findfirst(Base.Fix1(same_face,n2o_node[new_face]), old_f2n)
    @assert !isnothing(old)
    n2o_face[new] = old
  end

  return n2o_face
end

function _get_face_reindex(
  new_topo::PolytopalGridTopology, old_topo::UnstructuredGridTopology, d::Integer
)
  D = num_cell_dims(old_topo)
  cell_to_ctype = get_cell_type(old_topo)
  ctype_to_poly = get_polytopes(old_topo)
  ctype_to_faceperm = map(ctype_to_poly) do poly
    pos = _polytope_face_reindex(poly, d, true)
    neg = _polytope_face_reindex(poly, d, false)
    return [neg, pos]
  end

  cell_to_coords = get_cell_coordinates(old_topo)
  cell_polys = expand_cell_data(ctype_to_poly, cell_to_ctype)
  cell_to_pos = lazy_map(_polytope_is_positive, cell_polys, cell_to_coords)

  old_c2f = get_faces(old_topo, D, d)
  new_c2f = get_faces(new_topo, D, d)

  new_to_old_face = zeros(Int32, num_faces(new_topo, d))
  for (cell, ctype) in enumerate(cell_to_ctype)
    pos = cell_to_pos[cell]
    facemap = ctype_to_faceperm[ctype][pos+1]
    new_faces = dataview(new_c2f, cell)
    old_faces = dataview(old_c2f, cell)
    for (inew, iold) in enumerate(facemap)
      new = new_faces[inew]
      old = old_faces[iold]
      @check iszero(new_to_old_face[new]) || isequal(new_to_old_face[new], old)
      new_to_old_face[new] = old
    end
  end
  @check !any(iszero, new_to_old_face)

  return new_to_old_face
end

function generate_polytopes(
  cell_polys::AbstractArray{<:Polytope{2}},cell_nodes,node_coordinates
)
  o_cell_nodes = oriented_cell_nodes(cell_polys,cell_nodes,node_coordinates)
  cell_coordinates = lazy_map(Broadcasting(Reindex(node_coordinates)),o_cell_nodes)
  polytopes = map(Polygon,cell_coordinates)
  return polytopes, o_cell_nodes
end

function generate_polytopes(
  cell_polys::AbstractArray{<:Polytope{3}},cell_nodes,node_coordinates
)
  cell_coordinates = lazy_map(Broadcasting(Reindex(node_coordinates)),cell_nodes)
  polytopes = map(Polyhedron,cell_polys,cell_coordinates)
  return polytopes, cell_nodes
end

function compute_reffaces(::Type{Polytope{0}}, g::PolytopalGridTopology{Dc}) where Dc
  return [VERTEX], fill(Int8(1),num_faces(g,0))
end

function compute_reffaces(::Type{Polytope{1}}, g::PolytopalGridTopology{Dc}) where Dc
  [SEGMENT], fill(Int8(1),num_faces(g,1))
end

function compute_reffaces(::Type{Polytope{Dc}}, g::PolytopalGridTopology{Dc}) where Dc
  (get_polytopes(g), get_cell_type(g))
end

function compute_reffaces(::Type{Polytope{Df}}, g::PolytopalGridTopology{Dc}) where {Df,Dc}
  cell_to_face = get_faces(g,Dc,Df)
  face_to_cell = lazy_map(getindex, get_faces(g,Df,Dc), Fill(1,num_faces(g,Df)))
  face_to_lface = find_local_index(face_to_cell,cell_to_face)
  face_to_poly  = lazy_map(Reindex(get_polytopes(g)),face_to_cell)
  return map(Polytope{Df},face_to_poly,face_to_lface), Base.OneTo(num_faces(g,Df))
end

function num_faces(g::PolytopalGridTopology{Dc},d::Integer) where Dc
  if d == 0
    length(g.vertex_coordinates)
  elseif d == Dc
    length(g.polytopes)
  else
    face_to_cells = get_faces(g,d,Dc)
    length(face_to_cells)
  end
end

function get_faces(g::PolytopalGridTopology, dimfrom::Integer, dimto::Integer)
  _setup_faces!(g,dimfrom,dimto)
  g.n_m_to_nface_to_mfaces[dimfrom+1,dimto+1]
end

OrientationStyle(::Type{PolytopalGridTopology}) = NonOriented()
get_vertex_coordinates(g::PolytopalGridTopology) = g.vertex_coordinates
get_cell_type(g::PolytopalGridTopology) = Base.OneTo(length(g.polytopes))
get_polytopes(g::PolytopalGridTopology) = g.polytopes

function GridTopology(::Type{<:Polytope{Dc}},topo::PolytopalGridTopology{Dc}) where Dc
  return topo
end

function GridTopology(::Type{<:Polytope{Df}},topo::PolytopalGridTopology{Dc}) where {Df,Dc}
  node_coordinates = get_vertex_coordinates(topo)
  face_to_nodes = get_faces(topo,Df,0)
  if iszero(Df) || isone(Df)
    polytopes = [ifelse(iszero(Df), VERTEX, SEGMENT)]
    face_to_ftype = fill(Int8(1),num_faces(topo,Df))
    ftopo = UnstructuredGridTopology(node_coordinates, face_to_nodes, face_to_ftype, polytopes)
  else
    face_to_polytopes = get_reffaces(Polytope{Df},topo)
    ftopo = PolytopalGridTopology(node_coordinates, face_to_nodes, face_to_polytopes)
  end
  return ftopo
end

function restrict(topo::PolytopalGridTopology, cell_to_parent_cell::AbstractVector{<:Integer})
  D = num_cell_dims(topo)
  n_m_to_parent_nface_to_parent_mface = [
    Table(get_faces(topo,n,m)) for n in 0:D, m in 0:D 
  ]
  d_to_dface_to_parent_dface = [
    _setup_dface_to_parent_dface(
      Val{d}(),
      Val{D}(),
      n_m_to_parent_nface_to_parent_mface[D+1,d+1],
      num_faces(topo,d),
      cell_to_parent_cell
    ) for d in 0:D
  ]
  d_to_dface_to_vertices  = [
    _setup_connectivities_d(
      n_m_to_parent_nface_to_parent_mface[d+1,0+1],
      d_to_dface_to_parent_dface[d+1],
      d_to_dface_to_parent_dface[0+1],
      num_faces(topo,0)
    ) for d in 0:D
  ]
  vertex_coordinates = get_vertex_coordinates(topo)[d_to_dface_to_parent_dface[0+1]]
  polytopes = get_polytopes(topo)[cell_to_parent_cell]
  subtopo = PolytopalGridTopology(vertex_coordinates,d_to_dface_to_vertices,polytopes)
  return subtopo, d_to_dface_to_parent_dface
end

############################################################################################

"""
    struct PolytopalGrid{Dc,Dp,Tp,Tn} <: Grid{Dc,Dp}
      node_coordinates::Vector{Point{Dp,Tp}}
      cell_node_ids::Table{Int32,Vector{Int32},Vector{Int32}}
      polytopes::Vector{GeneralPolytope{Dc,Dp,Tp,Nothing}}
      facet_normal::Tn
    end

Grid for polytopal meshes.

Constructors: 

    PolytopalGrid(grid::Grid)
    PolytopalGrid(topo::PolytopalGridTopology)
    PolytopalGrid(node_coordinates,cell_node_ids,polytopes[,facet_normal=nothing])
    PolytopalGrid(polytopes::AbstractVector{<:GeneralPolytope})

"""
struct PolytopalGrid{Dc,Dp,Tp,Tn} <: Grid{Dc,Dp}
  node_coordinates::Vector{Point{Dp,Tp}}
  cell_node_ids::Table{Int32,Vector{Int32},Vector{Int32}}
  polytopes::Vector{<:GeneralPolytope{Dc,Dp,Tp,Nothing}}
  facet_normal::Tn

  function PolytopalGrid(
    node_coordinates::Vector{Point{Dp,Tp}},
    cell_node_ids::Table{Ti},
    polytopes::Vector{<:GeneralPolytope{Dc}},
    facet_normal=nothing;
  ) where {Dc,Dp,Tp,Ti}
    Tn = typeof(facet_normal)
    new{Dc,Dp,Tp,Tn}(
      node_coordinates, cell_node_ids, polytopes, facet_normal,
    )
  end
end

function PolytopalGrid(grid::Grid{Dc}) where Dc
  node_coordinates = collect1d(get_node_coordinates(grid))
  cell_polys = expand_cell_data(get_polytopes(grid),get_cell_type(grid))
  polytopes, cell_node_ids = generate_polytopes(
    cell_polys,Table(get_cell_node_ids(grid)),node_coordinates
  )
  PolytopalGrid(node_coordinates,cell_node_ids,polytopes)
end

function PolytopalGrid(topo::PolytopalGridTopology{Dc}) where Dc
  PolytopalGrid(get_vertex_coordinates(topo), get_faces(topo,Dc,0), get_polytopes(topo))
end

function PolytopalGrid(polytopes::AbstractVector{<:GeneralPolytope})
  node_coords = vcat((get_vertex_coordinates(q) for q in polytopes)...)
  offsets = [1,(cumsum(map(num_vertices,polytopes)).+1)...]
  cell_to_nodes = Table([ offsets[i]:offsets[i+1]-1 for i in 1:length(polytopes) ])
  return PolytopalGrid(node_coords, cell_to_nodes, polytopes)
end

is_first_order(::PolytopalGrid) = true
OrientationStyle(::Type{<:PolytopalGrid}) = NonOriented()
get_reffes(g::PolytopalGrid) = @notimplemented
get_polytopes(g::PolytopalGrid) = g.polytopes
get_cell_type(g::PolytopalGrid) = Base.OneTo(length(g.polytopes))
get_node_coordinates(g::PolytopalGrid) = g.node_coordinates
get_cell_node_ids(g::PolytopalGrid) = g.cell_node_ids
get_cell_map(g::PolytopalGrid) = Fill(GenericField(identity),num_cells(g))

function get_facet_normal(g::PolytopalGrid)
  @assert g.facet_normal != nothing "This Grid does not have information about normals."
  g.facet_normal
end

get_cell_ref_coordinates(g::PolytopalGrid) = get_cell_coordinates(g)

function restrict(grid::PolytopalGrid, cell_to_parent_cell::AbstractVector{Bool})
  cell_to_parent_cell = findall(cell_to_parent_cell)
  restrict(grid,cell_to_parent_cell)
end

function restrict(grid::PolytopalGrid, cell_to_parent_cell::AbstractVector{<:Integer})
  parent_cell_to_parent_nodes = get_cell_node_ids(grid)
  nparent_nodes = num_nodes(grid)
  node_to_parent_node, parent_node_to_node = _find_active_nodes(
    parent_cell_to_parent_nodes,cell_to_parent_cell,nparent_nodes
  )
  cell_to_nodes = _renumber_cell_nodes(
    parent_cell_to_parent_nodes,parent_node_to_node,cell_to_parent_cell
  )
  node_coordinates = get_node_coordinates(grid)[node_to_parent_node]
  polytopes = get_polytopes(grid)[cell_to_parent_cell]
  return PolytopalGrid(node_coordinates,cell_to_nodes,polytopes)
end

# TODO: This is not exactly the correct, but to do it correctly we have to deal 
# with  type dispatching correctly...
Base.view(a::Geometry.PolytopalGrid,b::AbstractArray) = Geometry.restrict(a,b)

function Grid(::Type{ReferenceFE{Df}},p::GeneralPolytope{Dc}) where {Df,Dc}
  vertex_coordinates = get_vertex_coordinates(p)
  face_to_vertices = Table(get_faces(p,Df,0))
  if iszero(Df) || isone(Df)
    reffes = [ifelse(iszero(Df), VERTEX1, SEG2)]
    face_to_ftype = fill(Int8(1),num_faces(p,Df))
    grid = UnstructuredGrid(vertex_coordinates, face_to_vertices, reffes, face_to_ftype)
  else
    face_to_polytopes = get_reffaces(Polytope{Df},p)
    grid = PolytopalGrid(vertex_coordinates, face_to_vertices, face_to_polytopes)
  end
  return grid
end

############################################################################################

"""
    struct PolytopalDiscreteModel{Dc,Dp,Tp,Tn} <: DiscreteModel{Dc,Dp}
      grid::PolytopalGrid{Dc,Dp,Tp,Tn}
      grid_topology::PolytopalGridTopology{Dc,Dp,Tp}
      labels::FaceLabeling
    end

Discrete model for polytopal grids.

Constructors: 

    PolytopalDiscreteModel(model::DiscreteModel)
    PolytopalDiscreteModel(grid::PolytopalGrid,grid_topology::PolytopalGridTopology,labels::FaceLabeling)
"""
struct PolytopalDiscreteModel{Dc,Dp,Tp,Tn} <: DiscreteModel{Dc,Dp}
  grid::PolytopalGrid{Dc,Dp,Tp,Tn}
  grid_topology::PolytopalGridTopology{Dc,Dp,Tp}
  labels::FaceLabeling
end

function PolytopalDiscreteModel(model::DiscreteModel)
  D = num_cell_dims(model)
  grid = PolytopalGrid(get_grid(model))
  topo = PolytopalGridTopology(get_grid_topology(model))
  d_to_dface_to_parent_dface = [
    _get_face_reindex(topo,get_grid_topology(model),d) for d in 0:D
  ]
  labels = restrict(get_face_labeling(model),d_to_dface_to_parent_dface)
  return PolytopalDiscreteModel(grid,topo,labels)
end

PolytopalDiscreteModel(model::PolytopalDiscreteModel) = model

get_grid(model::PolytopalDiscreteModel) = model.grid
get_grid_topology(model::PolytopalDiscreteModel) = model.grid_topology
get_face_labeling(model::PolytopalDiscreteModel) = model.labels

get_cell_map(model::PolytopalDiscreteModel) = get_cell_map(get_grid(model))

function compute_face_nodes(model::PolytopalDiscreteModel,d::Integer)
  topo = get_grid_topology(model)
  get_faces(topo,d,0)
end

function compute_face_own_nodes(model::PolytopalDiscreteModel,d::Integer)
  topo = get_grid_topology(model)
  if iszero(d)
    get_faces(topo,d,0)
  else
    empty_table(num_faces(topo,d))
  end
end

function Grid(::Type{ReferenceFE{Df}},model::PolytopalDiscreteModel{Dc}) where {Df,Dc}
  topo = get_grid_topology(model)
  node_coordinates = get_vertex_coordinates(topo)
  face_to_nodes = Table(get_faces(topo,Df,0))
  if iszero(Df) || isone(Df)
    reffes = [ifelse(iszero(Df), VERTEX1, SEG2)]
    face_to_ftype = fill(Int8(1),num_faces(topo,Df))
    grid = UnstructuredGrid(node_coordinates, face_to_nodes, reffes, face_to_ftype)
  else
    face_to_polytopes = get_reffaces(Polytope{Df},get_grid_topology(model))
    grid = PolytopalGrid(node_coordinates, face_to_nodes, face_to_polytopes)
  end
  return grid
end

function Grid(::Type{ReferenceFE{Dc}},model::PolytopalDiscreteModel{Dc}) where {Dc}
  get_grid(model)
end

function restrict(model::PolytopalDiscreteModel, cell_to_parent_cell::AbstractVector{<:Integer})
  grid = restrict(get_grid(model),cell_to_parent_cell)
  topo, d_to_dface_to_parent_dface = restrict(get_grid_topology(model),cell_to_parent_cell)
  labels = restrict(get_face_labeling(model),d_to_dface_to_parent_dface)
  return PolytopalDiscreteModel(grid,topo,labels)
end

@inline restrict(model::PolytopalDiscreteModel, ::IdentityVector) = model

############################################################################################
# Voronoi meshes

function orient_nodes_positive(v::AbstractVector{<:Point};keepfirst=false)
  _angle(c,v) = atan(v[1]-c[1],v[2]-c[2])
  perm = sortperm(v, by = Base.Fix1(_angle, mean(v)), rev = true)
  if keepfirst
    k = findfirst(==(1),perm)
    circshift!(perm,1 - k)
  end
  return perm
end

function is_positive_node_orientation(v::AbstractVector{<:Point})
  nv = length(v)
  c = mean(v)
  _angle(a,b) = atan(a[1]*b[2] - a[2]*b[1], dot(a,b))
  for i in eachindex(v)
    v1, v2 = v[i], v[mod(i,nv) + 1]
    if _angle(v1 - c, v2 - c) < 0
      return false
    end
  end
  return true
end

"""
    voronoi(model::DiscreteModel) -> PolytopalDiscreteModel
    voronoi(topo::GridTopology) -> PolytopalGridTopology

Given a mesh, computes it's associated Voronoi mesh. 
NOTE: Only working in 2D.
"""
function voronoi(model::DiscreteModel)
  new_topo = voronoi(get_grid_topology(model))
  new_grid = PolytopalGrid(new_topo)
  new_labels = FaceLabeling(new_topo) # TODO: Map labels from the original model
  return PolytopalDiscreteModel(new_grid,new_topo,new_labels)
end

function voronoi(topo::GridTopology{Dc}) where Dc
  n2c = get_faces(topo,0,Dc)
  n2e = get_faces(topo,0,Dc-1)
  c2n = get_faces(topo,Dc,0)
  e2n = get_faces(topo,Dc-1,0)
  node_coords = get_vertex_coordinates(topo)

  is_bnode = Geometry.get_isboundary_face(topo,0)
  is_bface = Geometry.get_isboundary_face(topo,Dc-1)
  bnodes = findall(is_bnode)
  bfaces = findall(is_bface)
  node_reindex = find_inverse_index_map(bnodes,num_faces(topo,0))    # Not counting the interior nodes for the global dual mesh
  face_reindex = find_inverse_index_map(bfaces,num_faces(topo,Dc-1))  # Not counting the interior faces for the global dual mesh    

  bulk_cen = map(mean,lazy_map(Broadcasting(Reindex(node_coords)), c2n))
  edge_cen = map(mean,lazy_map(Broadcasting(Reindex(node_coords)), e2n))

  new_node_coords = vcat(bulk_cen, edge_cen[bfaces], node_coords[bnodes])
  n_new_cells = num_faces(topo,0)

  ptrs = zeros(Int,n_new_cells+1)
  data = zeros(Int,n_new_cells * 3*2^(Dc-1)) # rough upper-bound

  f_offset = num_faces(topo,Dc)
  n_offset = f_offset + length(bfaces)

  ptrs[1] = 1
  for n in 1:n_new_cells
    if !is_bnode[n]
      # Case 1: Interior node
      #  --> New cell vertices are the centroids of the neighboring cells
      nbr_cells = view(n2c,n)
      ptrs[n+1] = ptrs[n] + length(nbr_cells)
      connectivity = nbr_cells
    else
      # Case 2: Boundary node
      #  --> New cell vertices given by: 
      #      1. The centroids of the neighboring cells
      #      2. The centroids of the neighboring faces
      #      3. The original node
      nbr_cells = view(n2c,n)
      nbr_faces = face_reindex[filter(f -> is_bface[f], view(n2e,n))]
      new_node  = node_reindex[n]
      connectivity = vcat(nbr_cells, nbr_faces .+ f_offset, new_node + n_offset)
    end
    # Orient the vertices counter-clockwise
    new_coords = view(new_node_coords,connectivity)
    perm = orient_nodes_positive(new_coords)
    # Store the connectivity
    ptrs[n+1] = ptrs[n] + length(connectivity)
    data[ptrs[n]:ptrs[n+1]-1] .= connectivity[perm]
  end
  resize!(data,ptrs[end]-1)
  new_c2n = Arrays.Table(data,ptrs)

  # Create polytopes
  Tp = eltype(eltype(new_node_coords))
  polytopes = Vector{GeneralPolytope{Dc,Dc,Tp,Nothing}}(undef,n_new_cells)
  for cell in 1:n_new_cells
    vertices = view(new_c2n,cell)
    polytopes[cell] = Polygon(new_node_coords[vertices])
  end

  return PolytopalGridTopology(new_node_coords,new_c2n,polytopes)
end
