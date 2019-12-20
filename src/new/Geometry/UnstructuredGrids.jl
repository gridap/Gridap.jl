
"""
    struct UnstructuredGrid{Dc,Dp,Tp,Ti,O} <: Grid{Dc,Dp}
      node_coordinates::Vector{Point{Dp,Tp}}
      cell_nodes::Table{Ti,Int32}
      reffes::Vector{<:NodalReferenceFE{Dc}}
      cell_types::Vector{Int8}
    end
"""
struct UnstructuredGrid{Dc,Dp,Tp,Ti,O} <: Grid{Dc,Dp}
  node_coordinates::Vector{Point{Dp,Tp}}
  cell_nodes::Table{Ti,Int32}
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
    new{Dc,Dp,Tp,Ti,B}(node_coordinates,cell_nodes,reffes,cell_types)
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
  ::Type{UnstructuredGrid{Dc,Dp,Tp,Ti,B}}) where {Dc,Dp,Tp,Ti,B} = Val{B}()

get_reffes(g::UnstructuredGrid) = g.reffes

get_cell_type(g::UnstructuredGrid) = g.cell_types

get_node_coordinates(g::UnstructuredGrid) = g.node_coordinates

get_cell_nodes(g::UnstructuredGrid) = g.cell_nodes


# From ReferenceFE

"""
    UnstructuredGrid(reffe::NodalReferenceFE)

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

function UnstructuredGrid(::Type{<:ReferenceFE{D}},p::Polytope{D}) where D
  order = 1
  reffe = LagrangianRefFE(Float64,p,order)
  UnstructuredGrid(reffe)
end

"""
    UnstructuredGrid(::Type{<:ReferenceFE{d}},p::Polytope) where d
"""
function UnstructuredGrid(::Type{<:ReferenceFE{d}},p::Polytope) where d
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


## Low dim grids
#
#"""
#    UnstructuredGrid(::Type{<:ReferenceFE{d}},trian::Grid) where d
#"""
#function UnstructuredGrid(::Type{<:ReferenceFE{d}},trian::Grid) where d
#  model = UnstructuredDiscreteModel(trian)
#  UnstructuredGrid(NodalReferenceFE{d},model)
#end
#
#function generate_cell_to_faces(d, grid::UnstructuredGrid, cell_to_vertices, vertex_to_cells)
#  reffes = get_reffes(grid)
#  polytopes = map(get_polytope,reffes)
#  cell_type_to_lface_to_lvertices = map( (p)->get_faces(p,d,0), polytopes )
#  cell_to_cell_type = get_cell_type(grid)
#
#  generate_cell_to_faces(
#    cell_to_vertices,
#    cell_type_to_lface_to_lvertices,
#    cell_to_cell_type,
#    vertex_to_cells)
#
#end
#
#function generate_cell_to_vertices(grid::UnstructuredGrid)
#  if has_straight_faces(grid)
#    cell_to_vertices = get_cell_nodes(grid)
#    vertex_to_node = collect(1:num_nodes(grid))
#    node_to_vertex = vertex_to_node
#  else
#    cell_to_nodes = get_cell_nodes(grid)
#    cell_to_cell_type = get_cell_type(grid)
#    reffes = get_reffes(grid)
#    cell_type_to_lvertex_to_lnode = map(get_vertex_node, reffes)
#    cell_to_vertices, vertex_to_node, node_to_vertex = _generate_cell_to_vertices(
#      cell_to_nodes,cell_to_cell_type,cell_type_to_lvertex_to_lnode)
#  end
#  (cell_to_vertices, vertex_to_node, node_to_vertex)
#end

