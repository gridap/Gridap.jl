
"""
    struct UnstructuredGridTopology{Dc,Dp,T,O} <: GridTopology{Dc,Dp}
      # private fields
    end
"""
struct UnstructuredGridTopology{Dc,Dp,T,O} <: GridTopology{Dc,Dp}
  vertex_coordinates::Vector{Point{Dp,T}}
  n_m_to_nface_to_mfaces::Matrix{Table{Int32,Vector{Int32},Vector{Int32}}}
  cell_type::Vector{Int8}
  polytopes::Vector{Polytope{Dc}}
  orientation_style::O
end

# Constructors

"""
    UnstructuredGridTopology(
      vertex_coordinates::Vector{<:Point},
      cell_vertices::Table,
      cell_type::Vector{<:Integer},
      polytopes::Vector{<:Polytope},
      orientation_style::OrientationStyle=NonOriented())
"""
function UnstructuredGridTopology(
  vertex_coordinates::Vector{<:Point},
  cell_vertices::Table,
  cell_type::Vector{<:Integer},
  polytopes::Vector{<:Polytope},
  orientation_style::OrientationStyle=NonOriented())

  D = num_dims(first(polytopes))
  n = D+1
  n_m_to_nface_to_mfaces = Matrix{Table{Int32,Vector{Int32},Vector{Int32}}}(undef,n,n)
  n_m_to_nface_to_mfaces[D+1,0+1] = cell_vertices
  vertex_cells = generate_cells_around(cell_vertices,length(vertex_coordinates))
  n_m_to_nface_to_mfaces[0+1,D+1] = vertex_cells

  P = eltype(vertex_coordinates)
  Dp = length(P)
  T = eltype(P)
  O = typeof(orientation_style)

  UnstructuredGridTopology{D,Dp,T,O}(
    vertex_coordinates,
    n_m_to_nface_to_mfaces,
    cell_type,
    polytopes,orientation_style)

end

"""
    UnstructuredGridTopology(
      vertex_coordinates::Vector{<:Point},
      d_to_dface_vertices::Vector{<:Table},
      cell_type::Vector{<:Integer},
      polytopes::Vector{<:Polytope},
      orientation_style::OrientationStyle=NonOriented())
"""
function UnstructuredGridTopology(
  vertex_coordinates::Vector{<:Point},
  d_to_dface_vertices::Vector{<:Table},
  cell_type::Vector{<:Integer},
  polytopes::Vector{<:Polytope},
  orientation_style::OrientationStyle=NonOriented())

  D = num_dims(first(polytopes))
  n = D+1
  n_m_to_nface_to_mfaces = Matrix{Table{Int32,Vector{Int32},Vector{Int32}}}(undef,n,n)
  nvertices = length(vertex_coordinates)
  for d in 0:D
    dface_to_vertices = d_to_dface_vertices[d+1]
    n_m_to_nface_to_mfaces[d+1,0+1] = dface_to_vertices
    vertex_to_dfaces = generate_cells_around(dface_to_vertices,nvertices)
    n_m_to_nface_to_mfaces[0+1,d+1] = vertex_to_dfaces
  end

  P = eltype(vertex_coordinates)
  Dp = length(P)
  T = eltype(P)
  O = typeof(orientation_style)

  UnstructuredGridTopology{D,Dp,T,O}(
    vertex_coordinates,
    n_m_to_nface_to_mfaces,
    cell_type,
    polytopes,
    orientation_style)

end

"""
    UnstructuredGridTopology(topo::GridTopology)
"""
function UnstructuredGridTopology(topo::GridTopology)

  vertex_coordinates = collect1d(get_vertex_coordinates(topo))
  cell_type = collect1d(get_cell_type(topo))
  polytopes = get_polytopes(topo)
  D = num_cell_dims(topo)
  d_to_dface_vertices = [ Table(get_faces(topo,d,0)) for d in 0:D ]
  orientation = OrientationStyle(topo)

  UnstructuredGridTopology(
    vertex_coordinates,
    d_to_dface_vertices,
    cell_type,
    polytopes,
    orientation)

end

function UnstructuredGridTopology(topo::UnstructuredGridTopology)
  topo
end

"""
    UnstructuredGridTopology(grid::UnstructuredGrid)

    UnstructuredGridTopology(
      grid::UnstructuredGrid,
      cell_to_vertices::Table,
      vertex_to_node::AbstractVector)
"""
function UnstructuredGridTopology(grid::UnstructuredGrid)
  cell_to_vertices, vertex_to_node, = _generate_cell_to_vertices_from_grid(grid)
  _generate_grid_topology_from_grid(grid,cell_to_vertices,vertex_to_node)
end

function UnstructuredGridTopology(grid::UnstructuredGrid, cell_to_vertices::Table, vertex_to_node::AbstractVector)
  _generate_grid_topology_from_grid(grid,cell_to_vertices,vertex_to_node)
end

function _generate_grid_topology_from_grid(grid::UnstructuredGrid,cell_to_vertices,vertex_to_node)

  @notimplementedif (! is_regular(grid)) "Extrtacting the GridTopology form a Grid only implemented for the regular case"

  node_to_coords = get_node_coordinates(grid)
  if vertex_to_node == 1:num_nodes(grid)
    vertex_to_coords = node_to_coords
  else
    vertex_to_coords = node_to_coords[vertex_to_node]
  end

  cell_to_type = get_cell_type(grid)
  polytopes = map(get_polytope, get_reffes(grid))

  UnstructuredGridTopology(
    vertex_to_coords,
    cell_to_vertices,
    cell_to_type,
    polytopes,
    OrientationStyle(grid))

end

function _generate_cell_to_vertices_from_grid(grid::UnstructuredGrid)
  if is_first_order(grid)
    cell_to_vertices = Table(get_cell_nodes(grid))
    vertex_to_node = collect(1:num_nodes(grid))
    node_to_vertex = vertex_to_node
  else
    cell_to_nodes = get_cell_nodes(grid)
    cell_to_cell_type = get_cell_type(grid)
    reffes = get_reffes(grid)
    cell_type_to_lvertex_to_lnode = map(get_vertex_node, reffes)
    cell_to_vertices, vertex_to_node, node_to_vertex = _generate_cell_to_vertices(
      cell_to_nodes,
      cell_to_cell_type,
      cell_type_to_lvertex_to_lnode,
      num_nodes(grid))
  end
  (cell_to_vertices, vertex_to_node, node_to_vertex)
end


function _generate_cell_to_vertices(
  cell_to_nodes::Table,
  cell_to_cell_type::AbstractVector{<:Integer},
  cell_type_to_lvertex_to_lnode::Vector{Vector{Int}},
  nnodes::Int=maximum(cell_to_nodes.data))

  data, ptrs, vertex_to_node, node_to_vertex = _generate_cell_to_vertices(
    cell_to_nodes.data,
    cell_to_nodes.ptrs,
    cell_to_cell_type,
    cell_type_to_lvertex_to_lnode,
    nnodes)

  (Table(data,ptrs), vertex_to_node)
end

function _generate_cell_to_vertices(
  cell_to_nodes_data,
  cell_to_nodes_ptrs,
  cell_to_cell_type,
  cell_type_to_lvertex_to_lnode,
  nnodes)

  cell_to_vertices_ptrs = similar(cell_to_nodes_ptrs)

  cell_type_to_nlvertices = map(length,cell_type_to_lvertex_to_lnode)

  _generate_cell_to_vertices_count!(
    cell_to_vertices_ptrs,
    cell_to_cell_type,
    cell_type_to_nlvertices)

  T = eltype(cell_to_nodes_data)

  node_to_vertex = fill(T(UNSET),nnodes)

  length_to_ptrs!(cell_to_vertices_ptrs)

  ndata = cell_to_vertices_ptrs[end]-1
  cell_to_vertices_data = zeros(T,ndata)

  _generate_cell_to_vertices_fill!(
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    cell_to_nodes_data,
    cell_to_nodes_ptrs,
    node_to_vertex,
    cell_to_cell_type,
    cell_type_to_lvertex_to_lnode)

  vertex_to_node = find_inverse_index_map(node_to_vertex)

  (cell_to_vertices_data, cell_to_vertices_ptrs, vertex_to_node, node_to_vertex)

end

function  _generate_cell_to_vertices_count!(
    cell_to_vertices_ptrs,
    cell_to_cell_type,
    cell_type_to_nlvertices)

  cells = 1:length(cell_to_cell_type)
  for cell in cells
    cell_type = cell_to_cell_type[cell]
    nlvertices = cell_type_to_nlvertices[cell_type]
    cell_to_vertices_ptrs[1+cell] = nlvertices
  end
end

function  _generate_cell_to_vertices_fill!(
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    cell_to_nodes_data,
    cell_to_nodes_ptrs,
    node_to_vertex,
    cell_to_cell_type,
    cell_type_to_lvertex_to_lnode)

  cells = 1:length(cell_to_cell_type)

  vertex = 1

  for cell in cells
    cell_type = cell_to_cell_type[cell]
    a = cell_to_nodes_ptrs[cell]-1
    b = cell_to_vertices_ptrs[cell]-1

    lvertex_to_lnode = cell_type_to_lvertex_to_lnode[cell_type]
    for (lvertex, lnode) in enumerate(lvertex_to_lnode)
      node = cell_to_nodes_data[a+lnode]
      if node_to_vertex[node] == UNSET
        node_to_vertex[node] = vertex
        vertex += 1
      end
      cell_to_vertices_data[b+lvertex] = node_to_vertex[node]
    end

  end
end

function GridTopology(::Type{<:Polytope{D}},topo::GridTopology{D}) where D
  topo
end

function GridTopology(::Type{<:Polytope{D}},topo::GridTopology) where D
  _topo = UnstructuredGridTopology(topo)
  GridTopology(Polytope{D},_topo)
end

function GridTopology(::Type{<:Polytope{D}},topo::UnstructuredGridTopology) where D

  vertex_coordinates = get_vertex_coordinates(topo)
  cell_type = get_face_type(topo,D)
  polytopes = get_reffaces(Polytope{D},topo)

  d_to_dface_vertices = [ get_faces(topo,d,0) for d in 0:D ]
  orientation = OrientationStyle(topo)

  topo_d = UnstructuredGridTopology(
    vertex_coordinates,
    d_to_dface_vertices,
    cell_type,
    polytopes,
    orientation)

end

# Needed, do not remove
function num_faces(g::UnstructuredGridTopology,d::Integer)
  if d == 0
    length(g.vertex_coordinates)
  elseif d == num_cell_dims(g)
    length(g.cell_type)
  else
    D = num_cell_dims(g)
    face_to_cells = get_faces(g,d,D)
    length(face_to_cells)
  end
end

# Implementation of abstract API

OrientationStyle(::Type{UnstructuredGridTopology{Dc,Dp,T,O}}) where {Dc,Dp,T,O} = O()

get_vertex_coordinates(g::UnstructuredGridTopology) = g.vertex_coordinates

get_cell_type(g::UnstructuredGridTopology) = g.cell_type

get_polytopes(g::UnstructuredGridTopology) = collect1d(g.polytopes)

function get_faces(g::UnstructuredGridTopology, dimfrom::Integer, dimto::Integer)
  _setup_faces!(g,dimfrom,dimto)
  g.n_m_to_nface_to_mfaces[dimfrom+1,dimto+1]
end

function _setup_faces!(g,dimfrom,dimto)
  if isassigned(g.n_m_to_nface_to_mfaces,dimfrom+1,dimto+1)
    return nothing
  end
  D = num_cell_dims(g)
  if dimfrom==0
    if dimto==0
      _setup_face_to_face!(g,dimfrom)
    elseif dimto==D
      nothing
    else
      _setup_face_to_vertices!(g,dimto)
    end
  elseif dimfrom == D
    if dimto==0
      nothing
    elseif dimto==D
      _setup_face_to_face!(g,dimfrom)
    else
      _setup_cell_to_faces!(g,dimto)
    end
  else
    if dimto==0
      _setup_face_to_vertices!(g,dimfrom)
    elseif dimto==D
      _setup_cell_to_faces!(g,dimfrom)
    elseif dimto==dimfrom
      _setup_face_to_face!(g,dimfrom)
    elseif dimfrom > dimto
      _setup_nface_to_mface!(g,dimfrom,dimto)
    else
      _setup_nface_to_mface!(g,dimto,dimfrom)
    end
  end
  nothing
end

function _setup_cell_to_faces!(model,dimto)

  D = num_cell_dims(model)

  if isassigned(model.n_m_to_nface_to_mfaces,D+1,dimto+1)
    return
  end

  if isassigned(model.n_m_to_nface_to_mfaces,dimto+1,0+1)

    cell_to_vertices = get_faces(model,D,0)
    vertex_to_faces = get_faces(model,0,dimto)
    cell_to_cell_type = get_cell_type(model)
    polytopes = get_polytopes(model)
    cell_type_to_lface_to_lvertices = map( (p)->get_faces(p,dimto,0), polytopes )

    cell_to_faces = find_cell_to_faces(
      cell_to_vertices,
      cell_type_to_lface_to_lvertices,
      cell_to_cell_type,
      vertex_to_faces)

  else

    cell_to_vertices = get_faces(model,D,0)
    vertex_to_cells = get_faces(model,0,D)
    cell_to_cell_type = get_cell_type(model)
    polytopes = get_polytopes(model)
    cell_type_to_lface_to_lvertices = map( (p)->get_faces(p,dimto,0), polytopes )

    cell_to_faces = generate_cell_to_faces(
        cell_to_vertices,
        cell_type_to_lface_to_lvertices,
        cell_to_cell_type,
        vertex_to_cells)

  end

  faces_to_cells = generate_cells_around(cell_to_faces)

  model.n_m_to_nface_to_mfaces[D+1,dimto+1] = cell_to_faces
  model.n_m_to_nface_to_mfaces[dimto+1,D+1] = faces_to_cells

  nothing

end

function _setup_face_to_vertices!(model,dimfrom)

  if isassigned(model.n_m_to_nface_to_mfaces,dimfrom+1,0+1)
    return
  end

  D = num_cell_dims(model)

  cell_to_vertices = get_faces(model,D,0)
  cell_to_faces = get_faces(model,D,dimfrom)
  cell_to_ctype = get_cell_type(model)
  polytopes = get_polytopes(model)
  ctype_to_lface_to_lvertices = map( (p)->get_faces(p,dimfrom,0), polytopes )

  nfaces = num_faces(model,dimfrom)

  face_to_vertices = generate_face_to_vertices(
    cell_to_vertices,
    cell_to_faces,
    cell_to_ctype,
    ctype_to_lface_to_lvertices,
    nfaces)

  nvertices = num_faces(model,0)

  vertex_to_faces = generate_cells_around(face_to_vertices,nvertices)

  model.n_m_to_nface_to_mfaces[dimfrom+1,0+1] = face_to_vertices
  model.n_m_to_nface_to_mfaces[0+1,dimfrom+1] = vertex_to_faces

  nothing

end

function _setup_face_to_face!(model,d)

  if isassigned(model.n_m_to_nface_to_mfaces,d+1,d+1)
    return
  end

  nfaces = num_faces(model,d)
  id = identity_table(Int32,Int32,nfaces)
  model.n_m_to_nface_to_mfaces[d+1,d+1] = id

  return

end

function _setup_nface_to_mface!(model,n,m)
  @assert n > m

  if isassigned(model.n_m_to_nface_to_mfaces,n+1,m+1)
    return
  end

  nface_to_vertices = get_faces(model,n,0)
  vertex_to_mfaces = get_faces(model,0,m)
  nface_to_nftype = get_face_type(model,n)

  polytopes = get_reffaces(Polytope{n},model)
  nftype_to_lmface_to_lvertices = map( (p)->get_faces(p,m,0), polytopes)

  nface_to_mfaces = find_cell_to_faces(
    nface_to_vertices,
    nftype_to_lmface_to_lvertices,
    nface_to_nftype,
    vertex_to_mfaces)

  num_mfaces = num_faces(model,m)
  mface_to_nfaces = generate_cells_around(nface_to_mfaces,num_mfaces)

  model.n_m_to_nface_to_mfaces[n+1,m+1] = nface_to_mfaces
  model.n_m_to_nface_to_mfaces[m+1,n+1] = mface_to_nfaces

  nothing

end
