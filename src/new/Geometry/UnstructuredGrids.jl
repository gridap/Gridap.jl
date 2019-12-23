
"""
    struct UnstructuredGrid{Dc,Dp,Tp,Ti,O} <: Grid{Dc,Dp}
      node_coordinates::Vector{Point{Dp,Tp}}
      cell_nodes::Table{Ti,Int32}
      reffes::Vector{<:NodalReferenceFE{Dc}}
      cell_types::Vector{Int8}
    end
"""
struct UnstructuredGrid{Dc,Dp,Tp,O} <: Grid{Dc,Dp}
  node_coordinates::Vector{Point{Dp,Tp}}
  cell_nodes::Table{Int,Int32}
  reffes::Vector{<:NodalReferenceFE{Dc}}
  cell_types::Vector{Int8}
  @doc """
      function UnstructuredGrid(
        node_coordinates::Vector{Point{Dp,Tp}},
        cell_nodes::Table{Ti},
        reffes::Vector{<:NodalReferenceFE{Dc}},
        cell_types::Vector,
        ::Val{B}=Val{false}()) where {Dc,Dp,Tp,Ti,B}
      end

  Low-level inner constructor.
  """
  function UnstructuredGrid(
    node_coordinates::Vector{Point{Dp,Tp}},
    cell_nodes::Table{Ti},
    reffes::Vector{<:NodalReferenceFE{Dc}},
    cell_types::Vector,
    ::Val{B}=Val{false}()) where {Dc,Dp,Tp,Ti,B}
    new{Dc,Dp,Tp,B}(node_coordinates,cell_nodes,reffes,cell_types)
  end
end

""" 
    UnstructuredGrid(trian::Grid)
"""
function UnstructuredGrid(trian::Grid)
  @assert is_regular(trian) "UnstructuredGrid constructor only for regular grids"
  node_coordinates = collect1d(get_node_coordinates(trian))
  cell_nodes = Table(get_cell_nodes(trian))
  reffes = get_reffes(trian)
  cell_types = collect1d(get_cell_type(trian))
  orien = OrientationStyle(trian)
  UnstructuredGrid(node_coordinates,cell_nodes,reffes,cell_types,orien)
end

function UnstructuredGrid(trian::UnstructuredGrid)
  trian
end

OrientationStyle(
  ::Type{UnstructuredGrid{Dc,Dp,Tp,B}}) where {Dc,Dp,Tp,B} = Val{B}()

get_reffes(g::UnstructuredGrid) = g.reffes

get_cell_type(g::UnstructuredGrid) = g.cell_types

get_node_coordinates(g::UnstructuredGrid) = g.node_coordinates

get_cell_nodes(g::UnstructuredGrid) = g.cell_nodes


# From ReferenceFE

"""
    UnstructuredGrid(reffe::LagrangianRefFE)

Build a grid with a single cell that is the given reference FE itself
"""
function UnstructuredGrid(reffe::LagrangianRefFE)
  node_coordinates = get_node_coordinates(reffe)
  cell_nodes = Table([collect(1:num_nodes(reffe)),])
  reffes = [reffe,]
  cell_types = [1,]
  UnstructuredGrid(node_coordinates,cell_nodes,reffes,cell_types)
end

# From Polytope

function UnstructuredGrid(::Type{ReferenceFE{D}},p::Polytope{D}) where D
  order = 1
  reffe = LagrangianRefFE(Float64,p,order)
  UnstructuredGrid(reffe)
end

"""
    UnstructuredGrid(::Type{ReferenceFE{d}},p::Polytope) where d
"""
function UnstructuredGrid(::Type{ReferenceFE{d}},p::Polytope) where d
  node_coordinates = get_vertex_coordinates(p)
  cell_nodes = Table(get_faces(p,d,0))
  reffaces = get_reffaces(Polytope{d},p)
  cell_type = get_face_type(p,d)
  order = 1
  reffes = map( (f)-> LagrangianRefFE(Float64,f,order), reffaces)
  UnstructuredGrid(
    node_coordinates,
    cell_nodes,
    reffes,
    cell_type)
end

# From coordinates

"""
    UnstructuredGrid(x::AbstractArray{<:Point})
"""
function UnstructuredGrid(x::AbstractArray{<:Point})
  np = length(x)
  node_coords = collect1d(x)
  cell_nodes = identity_table(Int,Int32,np)
  cell_type = fill(1,np)
  order = 1
  reffes = [LagrangianRefFE(Float64,VERTEX,order),]
  UnstructuredGrid(
    node_coords,
    cell_nodes,
    reffes,
    cell_type)
end

# Extract grid topology

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


