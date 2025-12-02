
"""
    struct UnstructuredGrid{Dc,Dp,Tp,Ti,O} <: Grid{Dc,Dp}
      node_coordinates::Vector{Point{Dp,Tp}}
      cell_node_ids::Table{Ti,Int32}
      reffes::Vector{<:LagrangianRefFE{Dc}}
      cell_types::Vector{Int8}
    end
"""
struct UnstructuredGrid{Dc,Dp,Tp,O,Tn} <: Grid{Dc,Dp}
  node_coordinates::Vector{Point{Dp,Tp}}
  cell_node_ids::Table{Int32,Vector{Int32},Vector{Int32}}
  reffes::Vector{LagrangianRefFE{Dc}}
  cell_types::Vector{Int8}
  orientation_style::O
  facet_normal::Tn
  cell_map
  @doc """
      function UnstructuredGrid(
        node_coordinates::Vector{Point{Dp,Tp}},
        cell_node_ids::Table{Ti},
        reffes::Vector{<:LagrangianRefFE{Dc}},
        cell_types::Vector,
        orientation_style::OrientationStyle=NonOriented()) where {Dc,Dp,Tp,Ti}
      end

  Low-level inner constructor.
  """
  function UnstructuredGrid(
    node_coordinates::Vector{Point{Dp,Tp}},
    cell_node_ids::Table{Ti},
    reffes::Vector{<:LagrangianRefFE{Dc}},
    cell_types::Vector,
    orientation_style::OrientationStyle=NonOriented(),
    facet_normal=nothing;
    has_affine_map=nothing) where {Dc,Dp,Tp,Ti}

    if has_affine_map === nothing
      _has_affine_map = get_has_affine_map(reffes)
    else
      _has_affine_map = has_affine_map
    end
    cell_map = _compute_cell_map(node_coordinates,cell_node_ids,reffes,cell_types,_has_affine_map)
    B = typeof(orientation_style)
    Tn = typeof(facet_normal)
    new{Dc,Dp,Tp,B,Tn}(
      node_coordinates,
      cell_node_ids,
      reffes,
      cell_types,
      orientation_style,
      facet_normal,
      cell_map)
  end
end

function get_has_affine_map(ctype_reffe)
  ctype_poly = map(get_polytope,ctype_reffe)
  has_affine_map =
    all(map(is_first_order,ctype_reffe)) &&
    ( all(map(is_simplex,ctype_poly)) || all(map(p->num_dims(p)==1,ctype_poly)) )
end

function _compute_cell_map(node_coords,cell_node_ids,ctype_reffe,cell_ctype, has_affine_map)
  cell_coords = lazy_map(Broadcasting(Reindex(node_coords)),cell_node_ids)
  ctype_shapefuns = map(get_shapefuns,ctype_reffe)
  cell_shapefuns = expand_cell_data(ctype_shapefuns,cell_ctype)
  default_cell_map = lazy_map(linear_combination,cell_coords,cell_shapefuns)
  ctype_poly = map(get_polytope,ctype_reffe)
  if has_affine_map
    ctype_q0 = map(p->zero(first(get_vertex_coordinates(p))),ctype_poly)
    cell_q0 = expand_cell_data(ctype_q0,cell_ctype)
    default_cell_grad = lazy_map(∇,default_cell_map)
    origins = lazy_map(evaluate,default_cell_map,cell_q0)
    gradients = lazy_map(evaluate,default_cell_grad,cell_q0)
    cell_map = lazy_map(Fields.affine_map,gradients,origins)
  else
    cell_map = default_cell_map
  end
  Fields.MemoArray(cell_map)
end

"""
    UnstructuredGrid(grid::Grid;kwargs...)
"""
function UnstructuredGrid(grid::Grid;kwargs...)
  @assert is_regular(grid) "UnstructuredGrid constructor only for regular grids"
  node_coordinates = collect1d(get_node_coordinates(grid))
  cell_node_ids = Table(get_cell_node_ids(grid))
  reffes = get_reffes(grid)
  cell_types = collect1d(get_cell_type(grid))
  orien = OrientationStyle(grid)
  UnstructuredGrid(node_coordinates,cell_node_ids,reffes,cell_types,orien;kwargs...)
end

function UnstructuredGrid(grid::UnstructuredGrid)
  grid
end

OrientationStyle(
  ::Type{<:UnstructuredGrid{Dc,Dp,Tp,B}}) where {Dc,Dp,Tp,B} = B()

get_reffes(g::UnstructuredGrid) = g.reffes

get_cell_type(g::UnstructuredGrid) = g.cell_types

get_node_coordinates(g::UnstructuredGrid) = g.node_coordinates

get_cell_node_ids(g::UnstructuredGrid) = g.cell_node_ids

get_cell_map(g::UnstructuredGrid) = g.cell_map

function get_facet_normal(g::UnstructuredGrid)
  @assert g.facet_normal != nothing "This Grid does not have information about normals."
  g.facet_normal
end

# From ReferenceFE

"""
    UnstructuredGrid(reffe::LagrangianRefFE)

Build a grid with a single cell that is the given reference FE itself
"""
function UnstructuredGrid(reffe::LagrangianRefFE)
  node_coordinates = get_node_coordinates(reffe)
  cell_node_ids = Table([collect(1:num_nodes(reffe)),])
  reffes = [reffe,]
  cell_types = [1,]
  UnstructuredGrid(node_coordinates,cell_node_ids,reffes,cell_types)
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
  cell_node_ids = Table(get_faces(p,d,0))
  reffaces = get_reffaces(Polytope{d},p)
  cell_type = get_face_type(p,d)
  order = 1
  reffes = map( (f)-> LagrangianRefFE(Float64,f,order), reffaces)
  UnstructuredGrid(
    node_coordinates,
    cell_node_ids,
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
  cell_node_ids = identity_table(Int32,Int32,np)
  cell_type = fill(1,np)
  order = 1
  reffes = [LagrangianRefFE(Float64,VERTEX,order),]
  UnstructuredGrid(
    node_coords,
    cell_node_ids,
    reffes,
    cell_type)
end

# IO

function to_dict(grid::UnstructuredGrid)
  dict = Dict{Symbol,Any}()
  x = get_node_coordinates(grid)
  dict[:node_coordinates] = reinterpret(eltype(eltype(x)),x)
  dict[:Dp] = num_point_dims(grid)
  dict[:cell_node_ids] = to_dict(get_cell_node_ids(grid))
  dict[:reffes] = map(to_dict, get_reffes(grid))
  dict[:cell_type] = get_cell_type(grid)
  dict[:orientation] = is_oriented(grid)
  dict
end

function from_dict(::Type{UnstructuredGrid},dict::Dict{Symbol,Any})
  x = dict[:node_coordinates]
  T = eltype(x)
  Dp = dict[:Dp]
  node_coordinates::Vector{Point{Dp,T}} = reinterpret(Point{Dp,T},x)
  cell_node_ids = from_dict(Table{Int32,Vector{Int32},Vector{Int32}},dict[:cell_node_ids])
  reffes = [ from_dict(LagrangianRefFE,reffe) for reffe in dict[:reffes]]
  cell_type::Vector{Int8} = dict[:cell_type]
  O::Bool = dict[:orientation]
  UnstructuredGrid(
    node_coordinates,
    cell_node_ids,
    reffes,
    cell_type,
    O ? Oriented() : NonOriented())
end
