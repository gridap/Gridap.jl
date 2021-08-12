"""
"""
struct ReferenceTriangulation{Dc,Dp,G} <: Triangulation{Dc,Dp}
  parent_trian::G
  function ReferenceTriangulation(t::Triangulation{Dc,Dp}) where {Dc,Dp}
     new{Dc,Dp,typeof(t)}(t)
  end
end

function have_compatible_domains(a::ReferenceTriangulation,b::Triangulation)
  have_compatible_domains(a.parent_trian,b)
end
function have_compatible_domains(a::Triangulation,b::ReferenceTriangulation)
  have_compatible_domains(a,b.parent_trian)
end
function have_compatible_domains(a::ReferenceTriangulation,b::ReferenceTriangulation)
  a === b
end

# Triangulation API
TriangulationStyle(::Type{<:ReferenceTriangulation}) = SubTriangulation()

get_background_triangulation(trian::ReferenceTriangulation) = get_background_triangulation(trian.parent_trian)

get_reffes(trian::ReferenceTriangulation) = get_reffes(trian.parent_trian)

function get_cell_coordinates(trian::ReferenceTriangulation)
  get_cell_ref_coordinates(trian.parent_trian)
end

function get_cell_type(trian::ReferenceTriangulation)
  get_cell_type(trian.parent_trian)
end

function get_cell_node_ids(trian::ReferenceTriangulation)
  get_cell_node_ids(trian.parent_trian)
end

function get_cell_map(trian::ReferenceTriangulation)
  cell_to_coords = get_cell_coordinates(trian)
  cell_to_shapefuns = get_cell_shapefuns(trian.parent_trian)
  lazy_map(linear_combination,cell_to_coords,cell_to_shapefuns)
end

function get_cell_to_bgcell(trian::ReferenceTriangulation)
  get_cell_to_bgcell(trian.parent_trian)
end

function get_facet_normal(trian::ReferenceTriangulation)
  get_facet_normal(trian.parent_trian)
end

function get_cell_ref_map(trian::ReferenceTriangulation)
  get_cell_ref_map(trian.parent_trian)
end
