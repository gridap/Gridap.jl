
"""
    struct UnstructuredGrid{Dc,Dp,Tp,Ti,O} <: Grid{Dc,Dp}
      node_coordinates::Vector{Point{Dp,Tp}}
      cell_nodes::Table{Ti,Int32}
      reffes::Vector{<:LagrangianRefFE{Dc}}
      cell_types::Vector{Int8}
    end
"""
struct UnstructuredGrid{Dc,Dp,Tp,O} <: Grid{Dc,Dp}
  node_coordinates::Vector{Point{Dp,Tp}}
  cell_nodes::Table{Int32,Vector{Int32},Vector{Int32}}
  reffes::Vector{LagrangianRefFE{Dc}}
  cell_types::Vector{Int8}
  orientation_style::O
  @doc """
      function UnstructuredGrid(
        node_coordinates::Vector{Point{Dp,Tp}},
        cell_nodes::Table{Ti},
        reffes::Vector{<:LagrangianRefFE{Dc}},
        cell_types::Vector,
        orientation_style::OrientationStyle=NonOriented()) where {Dc,Dp,Tp,Ti}
      end

  Low-level inner constructor.
  """
  function UnstructuredGrid(
    node_coordinates::Vector{Point{Dp,Tp}},
    cell_nodes::Table{Ti},
    reffes::Vector{<:LagrangianRefFE{Dc}},
    cell_types::Vector,
    orientation_style::OrientationStyle=NonOriented()) where {Dc,Dp,Tp,Ti}
    B = typeof(orientation_style)
    new{Dc,Dp,Tp,B}(node_coordinates,cell_nodes,reffes,cell_types,orientation_style)
  end
end

"""
    UnstructuredGrid(grid::Grid)
"""
function UnstructuredGrid(grid::Grid)
  @assert is_regular(grid) "UnstructuredGrid constructor only for regular grids"
  node_coordinates = collect1d(get_node_coordinates(grid))
  cell_nodes = Table(get_cell_nodes(grid))
  reffes = get_reffes(grid)
  cell_types = collect1d(get_cell_type(grid))
  orien = OrientationStyle(grid)
  UnstructuredGrid(node_coordinates,cell_nodes,reffes,cell_types,orien)
end

function UnstructuredGrid(grid::UnstructuredGrid)
  grid
end

OrientationStyle(
  ::Type{UnstructuredGrid{Dc,Dp,Tp,B}}) where {Dc,Dp,Tp,B} = B()

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
  cell_nodes = identity_table(Int32,Int32,np)
  cell_type = fill(1,np)
  order = 1
  reffes = [LagrangianRefFE(Float64,VERTEX,order),]
  UnstructuredGrid(
    node_coords,
    cell_nodes,
    reffes,
    cell_type)
end

# IO

function to_dict(grid::UnstructuredGrid)
  dict = Dict{Symbol,Any}()
  x = get_node_coordinates(grid)
  dict[:node_coordinates] = reinterpret(eltype(eltype(x)),x)
  dict[:Dp] = num_point_dims(grid)
  dict[:cell_nodes] = to_dict(get_cell_nodes(grid))
  dict[:reffes] = map(to_dict, get_reffes(grid))
  dict[:cell_type] = get_cell_type(grid)
  dict[:orientation] = is_oriented(grid)
  dict
end

function from_dict(::Type{UnstructuredGrid},dict::Dict{Symbol,Any})
  x = collect1d(dict[:node_coordinates])
  T = eltype(x)
  Dp = dict[:Dp]
  node_coordinates::Vector{Point{Dp,T}} = reinterpret(Point{Dp,T},x)
  cell_nodes = from_dict(Table{Int32,Vector{Int32},Vector{Int32}},dict[:cell_nodes])
  reffes = [ from_dict(LagrangianRefFE,reffe) for reffe in dict[:reffes]]
  cell_type::Vector{Int8} = dict[:cell_type]
  O::Bool = dict[:orientation]
  UnstructuredGrid(
    node_coordinates,
    cell_nodes,
    reffes,
    cell_type,
    O ? Oriented() : NonOriented())
end

function simplexify(grid::UnstructuredGrid)
  reffes = get_reffes(grid)
  @notimplementedif length(reffes) != 1
  reffe = first(reffes)
  order = 1
  @notimplementedif get_order(reffe) != order
  p = get_polytope(reffe)
  ltcell_to_lpoints, simplex = simplexify(p)
  cell_to_points = get_cell_nodes(grid)
  tcell_to_points = _refine_grid_connectivity(cell_to_points, ltcell_to_lpoints)
  ctype_to_reffe = [LagrangianRefFE(Float64,simplex,order),]
  tcell_to_ctype = fill(Int8(1),length(tcell_to_points))
  point_to_coords = get_node_coordinates(grid)
  UnstructuredGrid(
    point_to_coords,
    tcell_to_points,
    ctype_to_reffe,
    tcell_to_ctype,
    Oriented())
end

function _refine_grid_connectivity(cell_to_points::Table, ltcell_to_lpoints)
  data, ptrs = _refine_grid_connectivity(
    cell_to_points.data,
    cell_to_points.ptrs,
    ltcell_to_lpoints)
  Table(data,ptrs)
end

function _refine_grid_connectivity(
  cell_to_points_data::AbstractVector{T},
  cell_to_points_ptrs::AbstractVector{P},
  ltcell_to_lpoints) where {T,P}

  nltcells = length(ltcell_to_lpoints)
  ncells = length(cell_to_points_ptrs) - 1
  ntcells = ncells * nltcells

  tcell_to_points_ptrs = zeros(P,ntcells+1)

  _refine_grid_connectivity_count!(
    tcell_to_points_ptrs,
    ncells,
    ltcell_to_lpoints)

  length_to_ptrs!(tcell_to_points_ptrs)

  ndata = tcell_to_points_ptrs[end]-1

  tcell_to_points_data = zeros(T,ndata)

  _refine_grid_connectivity!(
    tcell_to_points_data,
    cell_to_points_data,
    cell_to_points_ptrs,
    ltcell_to_lpoints )

  (tcell_to_points_data, tcell_to_points_ptrs)

end

function  _refine_grid_connectivity_count!(
    tcell_to_points_ptrs,
    ncells,
    ltcell_to_lpoints)

  tcell = 1

  for cell in 1:ncells
    for lpoints in ltcell_to_lpoints
      tcell_to_points_ptrs[tcell+1] = length(lpoints)
      tcell +=1
    end
  end

end

function _refine_grid_connectivity!(
    tcell_to_points_data,
    cell_to_points_data,
    cell_to_points_ptrs,
    ltcell_to_lpoints )

  ncells = length(cell_to_points_ptrs) - 1

  k = 1
  for cell in 1:ncells
    a = cell_to_points_ptrs[cell]-1
    for lpoints in ltcell_to_lpoints
      for lpoint in lpoints
        point = cell_to_points_data[a+lpoint]
        tcell_to_points_data[k] = point
        k += 1
      end
    end
  end

end


