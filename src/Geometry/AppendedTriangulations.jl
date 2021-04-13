
function lazy_append(a::Triangulation,b::Triangulation)
  AppendedTriangulation(a,b)
end

struct AppendedTriangulation{Dc,Dp} <: Triangulation{Dc,Dp}
  a::Triangulation{Dc,Dp}
  b::Triangulation{Dc,Dp}
  function AppendedTriangulation(a::Triangulation{Dc,Dp},b::Triangulation{Dc,Dp}) where {Dc,Dp}
    @assert get_background_triangulation(a) === get_background_triangulation(b)
    new{Dc,Dp}(a,b)
  end
end

TriangulationStyle(::Type{<:AppendedTriangulation}) = SubTriangulation()

get_background_triangulation(trian::AppendedTriangulation) = get_background_triangulation(trian.a)

function get_cell_coordinates(trian::AppendedTriangulation)
  a = get_cell_coordinates(trian.a)
  b = get_cell_coordinates(trian.b)
  lazy_append(a,b)
end

function get_cell_ref_coordinates(trian::AppendedTriangulation)
  a = get_cell_ref_coordinates(trian.a)
  b = get_cell_ref_coordinates(trian.b)
  lazy_append(a,b)
end

function get_reffes(trian::AppendedTriangulation)
  vcat(get_reffes(trian.a),get_reffes(trian.b))
end

function get_cell_type(trian::AppendedTriangulation)
  a = get_cell_type(trian.a)
  b = get_cell_type(trian.b) .+ Int8(length(get_reffes(trian.a)))
  lazy_append(a,b)
end

#function reindex(f::AbstractArray,trian::AppendedTriangulation)
#  a = reindex(f,trian.a)
#  b = reindex(f,trian.b)
#  lazy_append(a,b)
#end

function get_cell_to_bgcell(trian::AppendedTriangulation)
  a = get_cell_to_bgcell(trian.a)
  b = get_cell_to_bgcell(trian.b)
  lazy_append(a,b)
end

function get_cell_ref_map(trian::AppendedTriangulation)
  a = get_cell_ref_map(trian.a)
  b = get_cell_ref_map(trian.b)
  lazy_append(a,b)
end

#function restrict(f::AbstractArray,trian::AppendedTriangulation)
#  a = restrict(f,trian.a)
#  b = restrict(f,trian.b)
#  lazy_append(a,b)
#end

# In this case, we do not want a lazy_append since it will become difficult to
# compress / expand the reffes.
#function get_cell_reffe(trian::AppendedTriangulation)
#  a = get_cell_reffe(trian.a)
#  b = get_cell_reffe(trian.b)
#  lazy_append(a,b)
#end

function get_cell_shapefuns(trian::AppendedTriangulation)
  a = get_cell_shapefuns(trian.a)
  b = get_cell_shapefuns(trian.b)
  lazy_append(a,b)
end

function get_cell_map(trian::AppendedTriangulation)
  a = get_cell_map(trian.a)
  b = get_cell_map(trian.b)
  lazy_append(a,b)
end

function get_facet_normal(trian::AppendedTriangulation)
  cm = get_cell_map(trian)
  a = get_facet_normal(trian.a)
  b = get_facet_normal(trian.b)
  lazy_append(a,b)
end

#function CellQuadrature(trian::AppendedTriangulation,degree1::Integer,degree2::Integer)
#  quad1 = CellQuadrature(trian.a,degree1)
#  quad2 = CellQuadrature(trian.b,degree2)
#  lazy_append(quad1,quad2)
#end
#
#function CellQuadrature(trian::AppendedTriangulation,degree::Integer)
#  CellQuadrature(trian,degree,degree)
#end

#function compress_contributions(cell_mat::AppendedArray,trian::AppendedTriangulation)
#  if length(cell_mat.a) == num_cells(trian.a) && length(cell_mat.b) == num_cells(trian.b)
#    a = compress_contributions(cell_mat.a,trian.a)
#    b = compress_contributions(cell_mat.b,trian.b)
#    return lazy_append(a,b)
#  else
#    return cell_mat
#  end
#end
#
#function compress_ids(cell_ids::AppendedArray,trian::AppendedTriangulation)
#  if length(cell_ids.a) == num_cells(trian.a) && length(cell_ids.b) == num_cells(trian.b)
#    a = compress_ids(cell_ids.a,trian.a)
#    b = compress_ids(cell_ids.b,trian.b)
#    return lazy_append(a,b)
#  else
#    return cell_ids
#  end
#end

function compress_contributions(cell_mat,trian::AppendedTriangulation)
  cell_to_bgcell = get_cell_to_bgcell(trian)
  ccell_mat = compress_contributions(cell_mat,cell_to_bgcell)
  ccell_mat
end

function compress_ids(cell_ids,trian::AppendedTriangulation)
  cell_to_bgcell = get_cell_to_bgcell(trian)
  compress_ids(cell_ids,cell_to_bgcell)
end
