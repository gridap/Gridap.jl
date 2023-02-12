
# Quadrature rule

"""
"""
struct CellQuadrature{DDS,IDS} <: CellDatum
  cell_quad::AbstractArray{<:Quadrature}
  cell_point::AbstractArray{<:AbstractArray{<:Point}}
  cell_weight::AbstractArray{<:AbstractArray{<:Real}}
  trian::Triangulation
  data_domain_style::DDS
  integration_domain_style::IDS
end

# Old constructor (for backward compatibility)
function CellQuadrature(
  cell_quad::AbstractArray{<:Quadrature},
  cell_point::AbstractArray{<:AbstractArray{<:Point}},
  cell_weight::AbstractArray{<:AbstractArray{<:Real}},
  trian::Triangulation,
  data_domain_style::DomainStyle)
  CellQuadrature(cell_quad,cell_point,cell_weight,trian,data_domain_style,PhysicalDomain())
end

function CellQuadrature(trian::Triangulation,quad::Tuple{<:QuadratureName,Any,Any})
  CellQuadrature(trian,quad,PhysicalDomain())
end


function CellQuadrature(trian::Triangulation,quad::Tuple{<:QuadratureName,Any,Any},ids::DomainStyle)
  name, args, kwargs = quad
  cell_quad = Quadrature(trian,name,args...;kwargs...)
  CellQuadrature(trian,cell_quad,ids)
end

function CellQuadrature(trian::Triangulation,degree::Integer;kwargs...)
  CellQuadrature(trian,degree,PhysicalDomain();kwargs...)
end

function CellQuadrature(trian::Triangulation,degree::Integer,ids::DomainStyle;kwargs...)
  cell_quad = Quadrature(trian,degree;kwargs...)
  CellQuadrature(trian,cell_quad,ids)
end

function CellQuadrature(trian::Triangulation,quad::Quadrature)
  CellQuadrature(trian,quad,PhysicalDomain())
end

function CellQuadrature(trian::Triangulation,quad::Quadrature,ids::DomainStyle)
  cell_quad = expand_cell_data([quad],get_cell_type(trian))
  CellQuadrature(trian,cell_quad,ids)
end

function CellQuadrature(trian::Triangulation,
  cell_quad::AbstractVector{<:Quadrature})
  CellQuadrature(trian,cell_quad,PhysicalDomain())
end

function CellQuadrature(trian::Triangulation,
                        cell_quad::AbstractVector{<:Quadrature},ids::DomainStyle)
  ctype_to_quad, cell_to_ctype = compress_cell_data(cell_quad)
  ctype_to_point = map(get_coordinates,ctype_to_quad)
  ctype_to_weight = map(get_weights,ctype_to_quad)
  cell_point = expand_cell_data(ctype_to_point,cell_to_ctype)
  cell_weight = expand_cell_data(ctype_to_weight,cell_to_ctype)
  CellQuadrature(cell_quad,cell_point,cell_weight,trian,ReferenceDomain(),ids)
end

function CellQuadrature(trian::AppendedTriangulation,degree1,degree2;kwargs...)
  CellQuadrature(trian,degree1,degree2,PhysicalDomain();kwargs...)
end

function CellQuadrature(trian::AppendedTriangulation,degree1,degree2,ids::DomainStyle;kwargs...)
  quad1 = CellQuadrature(trian.a,degree1,ids;kwargs...)
  quad2 = CellQuadrature(trian.b,degree2,ids;kwargs...)
  lazy_append(quad1,quad2,trian)
end

function CellQuadrature(trian::AppendedTriangulation,degree::Integer;kwargs...)
  CellQuadrature(trian,degree,degree,PhysicalDomain();kwargs...)
end

function CellQuadrature(trian::AppendedTriangulation,degree::Integer,ids::DomainStyle;kwargs...)
  CellQuadrature(trian,degree,degree,ids;kwargs...)
end

function CellQuadrature(trian::AppendedTriangulation,
  quad::Tuple{<:QuadratureName,Any,Any})
  CellQuadrature(trian,quad,quad,PhysicalDomain())
end

function CellQuadrature(trian::AppendedTriangulation,
                        quad::Tuple{<:QuadratureName,Any,Any},ids::DomainStyle)
  CellQuadrature(trian,quad,quad,ids)
end

function CellQuadrature(trian::AppendedTriangulation,quad::Quadrature)
  CellQuadrature(trian,quad,PhysicalDomain())
end

function CellQuadrature(trian::AppendedTriangulation,quad::Quadrature,ids::DomainStyle)
  CellQuadrature(trian,quad,quad,ids)
end

function lazy_append(
  quad1::CellQuadrature,
  quad2::CellQuadrature,
  trian::AppendedTriangulation=lazy_append(quad1.trian,quad2.trian))

  @notimplementedif DomainStyle(quad1) != DomainStyle(quad2)
  cell_quad = lazy_append(quad1.cell_quad,quad2.cell_quad)
  cell_point = lazy_append(quad1.cell_point,quad2.cell_point)
  cell_weight = lazy_append(quad1.cell_weight,quad2.cell_weight)
  CellQuadrature(cell_quad,cell_point,cell_weight,trian,DomainStyle(quad1),PhysicalDomain())
end

get_data(f::CellQuadrature) = f.cell_quad
get_triangulation(f::CellQuadrature) = f.trian
DomainStyle(::Type{CellQuadrature{DDS,IDS}}) where {DDS,IDS} = DDS()

function change_domain(a::CellQuadrature,::ReferenceDomain,::PhysicalDomain)
  @notimplemented
end

function change_domain(a::CellQuadrature,::PhysicalDomain,::ReferenceDomain)
  @notimplemented
end

function get_cell_points(a::CellQuadrature)
  CellPoint(a.cell_point,a.trian,a.data_domain_style)
end

function integrate(f::CellField,quad::CellQuadrature)
  trian_f = get_triangulation(f)
  trian_x = get_triangulation(quad)

  msg = """\n
    Your are trying to integrate a CellField using a CellQuadrature defined on incompatible
    triangulations. Verify that either the two objects are defined in the same triangulation
    or that the triangulaiton of the CellField is the background triangulation of the CellQuadrature.
    """
  @check is_change_possible(trian_f,trian_x) msg

  b = change_domain(f,quad.trian,quad.data_domain_style)
  x = get_cell_points(quad)
  bx = b(x)
  if quad.data_domain_style == PhysicalDomain() &&
            quad.integration_domain_style == PhysicalDomain()
    lazy_map(IntegrationMap(),bx,quad.cell_weight)
  elseif quad.data_domain_style == ReferenceDomain() &&
            quad.integration_domain_style == PhysicalDomain()
    cell_map = get_cell_map(quad.trian)
    cell_Jt = lazy_map(∇,cell_map)
    cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)
    lazy_map(IntegrationMap(),bx,quad.cell_weight,cell_Jtx)
  elseif quad.data_domain_style == ReferenceDomain() &&
            quad.integration_domain_style == ReferenceDomain()
    cell_map = Fill(GenericField(identity),length(bx))
    cell_Jt = lazy_map(∇,cell_map)
    cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)
    lazy_map(IntegrationMap(),bx,quad.cell_weight,cell_Jtx)
  else
    @notimplemented
  end
end

function integrate(a,quad::CellQuadrature)
  b = CellField(a,quad.trian,quad.data_domain_style)
  integrate(b,quad)
end

# Some syntactic sugar

struct Integrand
  object
end

const ∫ = Integrand

(*)(a::Integrand,b::CellQuadrature) = integrate(a.object,b)
(*)(b::CellQuadrature,a::Integrand) = integrate(a.object,b)

# Cell measure

function get_cell_measure(trian::Triangulation)
  quad = CellQuadrature(trian,0)
  cell_to_dV = integrate(1,quad)
end

function get_cell_measure(strian::Triangulation,ttrian::Triangulation)
  scell_measure = get_cell_measure(strian)
  move_contributions(scell_measure,strian,ttrian) |> collect
end

#"""
#Contributions added to the cells of the background Triangulation.
#"""
#function get_cell_measure(trian::Triangulation)
#  quad = CellQuadrature(trian,0)
#  cell_to_dV = integrate(1,quad)
#  model = get_background_model(trian)
#  bgtrian = Triangulation(model)
#  D = num_cell_dims(model)
#  glue = get_glue(trian,Val(D))
#  cell_to_bgcell = glue.tface_to_mface
#  bgcell_to_dV = zeros(num_cells(bgtrian))
#  _meas_K_fill!(bgcell_to_dV,cell_to_dV,cell_to_bgcell)
#  bgcell_to_dV
#end
#
#function _meas_K_fill!(bgcell_to_dV,cell_to_dV,cell_to_bgcell)
#  cache = array_cache(cell_to_dV)
#  for (cell, bgcell) in enumerate(cell_to_bgcell)
#    dV = getindex!(cache,cell_to_dV,cell)
#    bgcell_to_dV[bgcell] += dV
#  end
#end
