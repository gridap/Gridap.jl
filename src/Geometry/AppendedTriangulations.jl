
function lazy_append(a::Grid,b::Grid)
  AppendedGrid(a,b)
end

struct AppendedGrid{Dc,Dp,A,B} <: Grid{Dc,Dp}
  a::A
  b::B
  function AppendedGrid( a::Grid{Dc,Dp}, b::Grid{Dc,Dp}) where {Dc,Dp}
    new{Dc,Dp,typeof(a),typeof(b)}(a,b)
  end
end

function get_node_coordinates(t::AppendedGrid)
  xa = get_node_coordinates(t.a)
  xb = get_node_coordinates(t.b)
  if xa === xb
    xa
  else
    lazy_append(xa,xb)
  end
end

function get_cell_node_ids(t::AppendedGrid)
  xa = get_node_coordinates(t.a)
  xb = get_node_coordinates(t.b)
  idsa = get_cell_node_ids(t.a)
  idsb = get_cell_node_ids(t.b)
  if xa === xb
    lazy_append(idsa,idsb)
  else
    nxa = Int32(length(xa))
    idsc = lazy_map(Broadcasting(i->i+nxa),idsb)
    lazy_append(idsa,idsc)
  end
end

function get_reffes(t::AppendedGrid)
  ra = get_reffes(t.a)
  rb = get_reffes(t.b)
  if ra === rb
    ra
  else
    vcat(ra,rb)
  end
end

function get_cell_type(t::AppendedGrid)
  ra = get_reffes(t.a)
  rb = get_reffes(t.b)
  a = get_cell_type(t.a)
  b = get_cell_type(t.b)
  if ra === rb
    lazy_append(a,b)
  else
    nra = Int8(length(ra))
    c = lazy_map(i->i+nra,b)
    lazy_append(a,c)
  end
end

function get_cell_coordinates(trian::AppendedGrid)
  a = get_cell_coordinates(trian.a)
  b = get_cell_coordinates(trian.b)
  lazy_append(a,b)
end

function get_cell_ref_coordinates(trian::AppendedGrid)
  a = get_cell_ref_coordinates(trian.a)
  b = get_cell_ref_coordinates(trian.b)
  lazy_append(a,b)
end

# In this case, we do not want a lazy_append since it will become difficult to
# compress / expand the reffes.
#function get_cell_reffe(trian::AppendedGrid)
#  a = get_cell_reffe(trian.a)
#  b = get_cell_reffe(trian.b)
#  lazy_append(a,b)
#end

function get_cell_shapefuns(trian::AppendedGrid)
  a = get_cell_shapefuns(trian.a)
  b = get_cell_shapefuns(trian.b)
  lazy_append(a,b)
end

function get_cell_map(trian::AppendedGrid)
  a = get_cell_map(trian.a)
  b = get_cell_map(trian.b)
  lazy_append(a,b)
end

function get_facet_normal(trian::AppendedGrid)
  cm = get_cell_map(trian)
  a = get_facet_normal(trian.a)
  b = get_facet_normal(trian.b)
  lazy_append(a,b)
end

function lazy_append(a::Triangulation,b::Triangulation)
  AppendedTriangulation(a,b)
end

"""
    AppendedTriangulation(a::Triangulation{Dc,Dp}, b::Triangulation{Dc,Dp})

Union of two triangulations built on the same [`DiscreteModel`](@ref).
"""
struct AppendedTriangulation{Dc,Dp,A,B} <: Triangulation{Dc,Dp}
  a::A
  b::B
  function AppendedTriangulation(
    a::Triangulation{Dc,Dp}, b::Triangulation{Dc,Dp}) where {Dc,Dp}
    @assert get_background_model(a) === get_background_model(b)
    new{Dc,Dp,typeof(a),typeof(b)}(a,b)
  end
end

get_background_model(t::AppendedTriangulation) = get_background_model(t.a)

function get_grid(t::AppendedTriangulation)
  a = get_grid(t.a)
  b = get_grid(t.b)
  lazy_append(a,b)
end

function get_glue(t::AppendedTriangulation,::Val{D}) where D
  a = get_glue(t.a,Val(D))
  b = get_glue(t.b,Val(D))
  if a==nothing || b==nothing
    return nothing
  end
  lazy_append(a,b)
end

function lazy_append(a::FaceToFaceGlue,b::FaceToFaceGlue)
  tface_to_mface = lazy_append(a.tface_to_mface,b.tface_to_mface)
  tface_to_mface_map = lazy_append(a.tface_to_mface_map,b.tface_to_mface_map)
  mface_to_tface = nothing
  FaceToFaceGlue(tface_to_mface,tface_to_mface_map,mface_to_tface)
end

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

function move_contributions(
  cell_mat::AppendedArray, trian::AppendedTriangulation)

  if length(cell_mat.a) == num_cells(trian.a) && length(cell_mat.b) == num_cells(trian.b)
    a,ta = move_contributions(cell_mat.a,trian.a)
    b,tb = move_contributions(cell_mat.b,trian.b)
    return lazy_append(a,b), lazy_append(ta,tb)
  else
    return cell_mat, trian
  end
end

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
