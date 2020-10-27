
"""
"""
struct RestrictedTriangulation{Dc,Dp,G,A} <: Triangulation{Dc,Dp}
  parent_trian::G
  cell_to_parent_cell::A
  @doc """
  """
  function RestrictedTriangulation(
    parent_trian::Triangulation,
    cell_to_parent_cell::AbstractVector{<:Integer})

    Dc = num_cell_dims(parent_trian)
    Dp = num_point_dims(parent_trian)
    G = typeof(parent_trian)
    A = typeof(cell_to_parent_cell)
    new{Dc,Dp,G,A}(parent_trian,cell_to_parent_cell)
  end
end

function RestrictedTriangulation(
  parent_trian::Triangulation,
  parent_cell_to_mask::AbstractVector{Bool})

  cell_to_parent_cell = findall(parent_cell_to_mask)
  RestrictedTriangulation(parent_trian,cell_to_parent_cell)
end

TriangulationStyle(::Type{<:RestrictedTriangulation}) = SubTriangulation()

get_background_triangulation(trian::RestrictedTriangulation) = get_background_triangulation(trian.parent_trian)

get_node_coordinates(trian::RestrictedTriangulation) = get_node_coordinates(trian.parent_trian)

get_reffes(trian::RestrictedTriangulation) = get_reffes(trian.parent_trian)

function get_cell_coordinates(trian::RestrictedTriangulation)
  parent_cell_data = get_cell_coordinates(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_type(trian::RestrictedTriangulation)
  parent_cell_data = get_cell_type(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_nodes(trian::RestrictedTriangulation)
  parent_cell_data = get_cell_nodes(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_map(trian::RestrictedTriangulation)
  parent_cell_data = get_cell_map(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_id(trian::RestrictedTriangulation)
  parent_cell_data = get_cell_id(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_facet_normal(trian::RestrictedTriangulation)
  parent_cell_data = get_facet_normal(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_ref_map(trian::RestrictedTriangulation)
  parent_cell_data = get_cell_ref_map(trian.parent_trian) 
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

#function reindex(trian::Triangulation,indices)
#  RestrictedTriangulation(trian,collect(Int,indices))
#end
#
#function reindex(trian::Triangulation,indices::Vector{Int})
#  RestrictedTriangulation(trian,indices)
#end
#
#function reindex(trian::Triangulation,indices::IdentityVector)
#  trian
#end
#
#"""
#"""
#struct RestrictedTriangulation{Dc,Dp,G} <: Triangulation{Dc,Dp}
#  oldtrian::G
#  cell_to_oldcell::Vector{Int}
#  @doc """
#  """
#  function RestrictedTriangulation(oldtrian::Triangulation{Dc,Dp},cell_to_oldcell::Vector{Int}) where {Dc,Dp}
#    new{Dc,Dp,typeof(oldtrian)}(oldtrian,cell_to_oldcell)
#  end
#end
#
#function get_reffes(trian::RestrictedTriangulation)
#  get_reffes(trian.oldtrian)
#end
#
#function get_cell_type(trian::RestrictedTriangulation)
#  reindex(get_cell_type(trian.oldtrian),trian.cell_to_oldcell)
#end
#
#function get_cell_coordinates(trian::RestrictedTriangulation)
#  reindex(get_cell_coordinates(trian.oldtrian),trian.cell_to_oldcell)
#end
#
#function get_cell_map(trian::RestrictedTriangulation)
#  cell_map = get_cell_map(trian.oldtrian)
#  reindex(cell_map,trian.cell_to_oldcell)
#end
#
#function get_cell_id(trian::RestrictedTriangulation)
#  reindex(get_cell_id(trian.oldtrian),trian.cell_to_oldcell)
#end
#
#function restrict(f::AbstractArray,trian::RestrictedTriangulation)
#  reindex(restrict(f,trian.oldtrian),trian.cell_to_oldcell)
#end
#
#function get_normal_vector(trian::RestrictedTriangulation)
#  reindex(get_normal_vector(trian.oldtrian),trian.cell_to_oldcell)
#end
