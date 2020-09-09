
function lazy_append(a::Triangulation,b::Triangulation)
  AppendedTriangulation(a,b)
end

struct AppendedTriangulation{Dc,Dp} <: Triangulation{Dc,Dp}
  a::Triangulation{Dc,Dp}
  b::Triangulation{Dc,Dp}
  memo::Dict
  function AppendedTriangulation(a::Triangulation{Dc,Dp},b::Triangulation{Dc,Dp}) where {Dc,Dp}
    new{Dc,Dp}(a,b,Dict())
  end
end

get_memo(a::AppendedTriangulation) = a.memo

function get_cell_coordinates(trian::AppendedTriangulation)
  a = get_cell_coordinates(trian.a)
  b = get_cell_coordinates(trian.b)
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

function reindex(f::AbstractArray,trian::AppendedTriangulation)
  a = reindex(f,trian.a)
  b = reindex(f,trian.b)
  lazy_append(a,b)
end

function get_cell_id(trian::AppendedTriangulation)
  a = get_cell_id(trian.a)
  b = get_cell_id(trian.b)
  lazy_append(a,b)
end

function get_cell_reffes(trian::AppendedTriangulation)
  a = get_cell_reffes(trian.a)
  b = get_cell_reffes(trian.b)
  lazy_append(a,b)
end

function get_cell_shapefuns(trian::AppendedTriangulation)
  a = get_cell_shapefuns(trian.a)
  b = get_cell_shapefuns(trian.b)
  lazy_append(a,b)
end

function compute_cell_map(trian::AppendedTriangulation)
  a = get_cell_map(trian.a)
  b = get_cell_map(trian.b)
  lazy_append(a,b)
end

function get_normal_vector(trian::AppendedTriangulation)
  n1 = get_normal_vector(trian.a)
  n2 = get_normal_vector(trian.b)
  ϕ1 = get_cell_map(trian.a)
  ϕ2 = get_cell_map(trian.b)
  ϕ = get_cell_map(trian)
  lazy_append(n1∘ϕ1,n2∘ϕ2)∘inverse_map(ϕ)
end

function reference_cell_quadrature(trian::AppendedTriangulation,degree1::Integer,degree2::Integer)
  quad1 = reference_cell_quadrature(trian.a,degree1)
  quad2 = reference_cell_quadrature(trian.b,degree2)
  lazy_append(quad1,quad2)
end

function reference_cell_quadrature(trian::AppendedTriangulation,degree::Integer)
  reference_cell_quadrature(trian,degree,degree)
end

function CellQuadrature(trian::AppendedTriangulation,degree1::Integer,degree2::Integer)
  quad_ref = reference_cell_quadrature(trian,degree1,degree2)
  ϕ = get_cell_map(trian)
  ϕ(quad_ref)
end

function CellQuadrature(trian::AppendedTriangulation,degree::Integer)
  CellQuadrature(trian,degree,degree)
end

