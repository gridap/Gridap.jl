
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

# Constructors

function RestrictedTriangulation(
  parent_trian::Triangulation,
  parent_cell_to_mask::AbstractArray{Bool})

  cell_to_parent_cell = findall(collect1d(parent_cell_to_mask))
  RestrictedTriangulation(parent_trian,cell_to_parent_cell)
end

function RestrictedTriangulation(
  parent_trian::Triangulation,
  parent_cell_to_mask::AbstractVector{Bool})

  cell_to_parent_cell = findall(parent_cell_to_mask)
  RestrictedTriangulation(parent_trian,cell_to_parent_cell)
end

function Triangulation(
  parent_trian::Triangulation,
  cell_to_parent_cell::AbstractVector{<:Integer})

  RestrictedTriangulation(parent_trian,cell_to_parent_cell)
end

function Triangulation(
  parent_trian::Triangulation,
  parent_cell_to_mask::AbstractArray{Bool})

  RestrictedTriangulation(parent_trian,parent_cell_to_mask)
end

function Triangulation(
  parent_model::DiscreteModel,
  cell_to_parent_cell::AbstractVector{<:Integer})

  parent_trian = Triangulation(parent_model)
  RestrictedTriangulation(parent_trian,cell_to_parent_cell)
end

function Triangulation(
  parent_model::DiscreteModel,
  parent_cell_to_mask::AbstractArray{Bool})

  parent_trian = Triangulation(parent_model)
  RestrictedTriangulation(parent_trian,parent_cell_to_mask)
end

function Triangulation(parent_model::DiscreteModel, labels::FaceLabeling; tags)
  parent_trian = Triangulation(parent_model)
  parent_cell_to_mask = get_face_mask(labels,tags,num_cell_dims(parent_model))
  RestrictedTriangulation(parent_trian,parent_cell_to_mask)
end

# Triangulation API


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

