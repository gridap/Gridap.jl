
"""
"""
struct TriangulationPortion{Dc,Dp,G,A} <: Triangulation{Dc,Dp}
  parent_trian::G
  cell_to_parent_cell::A
  @doc """
  """
  function TriangulationPortion(
    parent_trian::Triangulation,
    cell_to_parent_cell::AbstractVector{<:Integer})

    Dc = num_cell_dims(parent_trian)
    Dp = num_point_dims(parent_trian)
    G = typeof(parent_trian)
    A = typeof(cell_to_parent_cell)
    new{Dc,Dp,G,A}(parent_trian,cell_to_parent_cell)
  end
end

TriangulationStyle(::Type{<:TriangulationPortion}) = SubTriangulation()

get_background_triangulation(trian::TriangulationPortion) = get_background_triangulation(trian.parent_trian)

get_node_coordinates(trian::TriangulationPortion) = get_node_coordinates(trian.parent_trian)

get_reffes(trian::TriangulationPortion) = get_reffes(trian.parent_trian)

function get_cell_coordinates(trian::TriangulationPortion)
  parent_cell_data = get_cell_coordinates(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_type(trian::TriangulationPortion)
  parent_cell_data = get_cell_type(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_nodes(trian::TriangulationPortion)
  parent_cell_data = get_cell_nodes(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_map(trian::TriangulationPortion)
  parent_cell_data = get_cell_map(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_id(trian::TriangulationPortion)
  parent_cell_data = get_cell_id(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_facet_normal(trian::TriangulationPortion)
  parent_cell_data = get_facet_normal(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_ref_map(trian::TriangulationPortion)
  parent_cell_data = get_cell_ref_map(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

#function reindex(trian::Triangulation,indices)
#  TriangulationPortion(trian,collect(Int,indices))
#end
#
#function reindex(trian::Triangulation,indices::Vector{Int})
#  TriangulationPortion(trian,indices)
#end
#
#function reindex(trian::Triangulation,indices::IdentityVector)
#  trian
#end
#
#"""
#"""
#struct TriangulationPortion{Dc,Dp,G} <: Triangulation{Dc,Dp}
#  oldtrian::G
#  cell_to_oldcell::Vector{Int}
#  @doc """
#  """
#  function TriangulationPortion(oldtrian::Triangulation{Dc,Dp},cell_to_oldcell::Vector{Int}) where {Dc,Dp}
#    new{Dc,Dp,typeof(oldtrian)}(oldtrian,cell_to_oldcell)
#  end
#end
#
#function get_reffes(trian::TriangulationPortion)
#  get_reffes(trian.oldtrian)
#end
#
#function get_cell_type(trian::TriangulationPortion)
#  reindex(get_cell_type(trian.oldtrian),trian.cell_to_oldcell)
#end
#
#function get_cell_coordinates(trian::TriangulationPortion)
#  reindex(get_cell_coordinates(trian.oldtrian),trian.cell_to_oldcell)
#end
#
#function get_cell_map(trian::TriangulationPortion)
#  cell_map = get_cell_map(trian.oldtrian)
#  reindex(cell_map,trian.cell_to_oldcell)
#end
#
#function get_cell_id(trian::TriangulationPortion)
#  reindex(get_cell_id(trian.oldtrian),trian.cell_to_oldcell)
#end
#
#function restrict(f::AbstractArray,trian::TriangulationPortion)
#  reindex(restrict(f,trian.oldtrian),trian.cell_to_oldcell)
#end
#
#function get_normal_vector(trian::TriangulationPortion)
#  reindex(get_normal_vector(trian.oldtrian),trian.cell_to_oldcell)
#end
