
"""
    struct UnstructuredGrid{Dc,Dp,Tp,Ti} <: ConformingTriangulation{Dc,Dp}
      node_coordinates::Vector{Point{Dp,Tp}}
      cell_nodes::Table{Ti,Int32}
      reffes::Vector{<:NodalReferenceFE{Dc}}
      cell_types::Vector{Int8}
    end
"""
struct UnstructuredGrid{Dc,Dp,Tp,Ti,B} <: ConformingTriangulation{Dc,Dp}
  node_coordinates::Vector{Point{Dp,Tp}}
  cell_nodes::Table{Ti,Int32}
  reffes::Vector{<:NodalReferenceFE{Dc}}
  cell_types::Vector{Int8}
  @doc """
      function UnstructuredGrid(
        node_coordinates::Vector{Point{Dp,Tp}},
        cell_nodes::Table{Ti},
        reffes::Vector{<:NodalReferenceFE{Dc}},
        cell_types::Vector) where {Dc,Dp,Tp,Ti}
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
    UnstructuredGrid(trian::ConformingTriangulation)
"""
function UnstructuredGrid(trian::ConformingTriangulation)
  @assert ConformityStyle(trian) == RegularConformity() "UnstructuredGrid constructor only for regular grids"
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
function UnstructuredGrid(reffe::NodalReferenceFE)
  node_coordinates = get_node_coordinates(reffe)
  cell_nodes = Table([collect(1:num_nodes(reffe)),])
  reffes = [reffe,]
  cell_types = [1,]
  UnstructuredGrid(node_coordinates,cell_nodes,reffes,cell_types)
end

# From Polytope

function UnstructuredGrid(::Type{<:ReferenceFE{D}},p::Polytope{D}) where D
  reffe = NodalReferenceFE(p)
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
  reffes = map(NodalReferenceFE,reffaces)
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
  reffes = [NodalReferenceFE(VERTEX),]
  UnstructuredGrid(
    node_coords,
    cell_nodes,
    reffes,
    cell_type)
end


# Low dim grids

"""
    UnstructuredGrid(::Type{<:ReferenceFE{d}},trian::ConformingTriangulation) where d
"""
function UnstructuredGrid(::Type{<:ReferenceFE{d}},trian::ConformingTriangulation) where d
  model = UnstructuredDiscreteModel(trian)
  UnstructuredGrid(NodalReferenceFE{d},model)
end

