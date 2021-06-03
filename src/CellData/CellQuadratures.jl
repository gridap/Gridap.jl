
# Quadrature rule

"""
"""
struct CellQuadrature{DS} <: CellDatum
  cell_quad::AbstractArray{<:Quadrature}
  cell_point::AbstractArray{<:AbstractArray{<:Point}}
  cell_weight::AbstractArray{<:AbstractArray{<:Real}}
  trian::Triangulation
  domain_style::DS
end


function CellQuadrature(trian::Triangulation,quad::Tuple{<:QuadratureName,Any,Any})
  name, args, kwargs = quad
  cell_quad = Quadrature(trian,name,args...;kwargs...)
  CellQuadrature(trian,cell_quad)
end

function CellQuadrature(trian::Triangulation,degree::Integer)
  cell_quad = Quadrature(trian,degree)
  CellQuadrature(trian,cell_quad)
end

function CellQuadrature(trian::Triangulation,quad::Quadrature)
  cell_quad = expand_cell_data([quad],get_cell_type(trian))
  CellQuadrature(trian,cell_quad)
end

function CellQuadrature(trian::Triangulation,cell_quad::AbstractVector{<:Quadrature})
  ctype_to_quad, cell_to_ctype = compress_cell_data(cell_quad)
  ctype_to_point = map(get_coordinates,ctype_to_quad)
  ctype_to_weight = map(get_weights,ctype_to_quad)
  cell_point = expand_cell_data(ctype_to_point,cell_to_ctype)
  cell_weight = expand_cell_data(ctype_to_weight,cell_to_ctype)
  CellQuadrature(cell_quad,cell_point,cell_weight,trian,ReferenceDomain())
end

function CellQuadrature(trian::AppendedTriangulation,degree1,degree2)
  quad1 = CellQuadrature(trian.a,degree1)
  quad2 = CellQuadrature(trian.b,degree2)
  lazy_append(quad1,quad2,trian)
end

function CellQuadrature(trian::AppendedTriangulation,degree::Integer)
  CellQuadrature(trian,degree,degree)
end

function CellQuadrature(trian::AppendedTriangulation,quad::Tuple{<:QuadratureName,Any,Any})
  CellQuadrature(trian,quad,quad)
end

function CellQuadrature(trian::AppendedTriangulation,quad::Quadrature)
  CellQuadrature(trian,quad,quad)
end

function lazy_append(
  quad1::CellQuadrature,
  quad2::CellQuadrature,
  trian::AppendedTriangulation=lazy_append(quad1.trian,quad2.trian))

  @notimplementedif DomainStyle(quad1) != DomainStyle(quad2)
  cell_quad = lazy_append(quad1.cell_quad,quad2.cell_quad)
  cell_point = lazy_append(quad1.cell_point,quad2.cell_point)
  cell_weight = lazy_append(quad1.cell_weight,quad2.cell_weight)
  CellQuadrature(cell_quad,cell_point,cell_weight,trian,DomainStyle(quad1))
end

get_data(f::CellQuadrature) = f.cell_quad
get_triangulation(f::CellQuadrature) = f.trian
DomainStyle(::Type{CellQuadrature{DS}}) where DS = DS()

function change_domain(a::CellQuadrature,::ReferenceDomain,::PhysicalDomain)
  @notimplemented
end

function change_domain(a::CellQuadrature,::PhysicalDomain,::ReferenceDomain)
  @notimplemented
end

function get_cell_points(a::CellQuadrature)
  CellPoint(a.cell_point,a.trian,a.domain_style)
end

function integrate(f::CellField,quad::CellQuadrature)

  trian_f = get_triangulation(f)
  trian_x = get_triangulation(quad)

  if have_compatible_domains(trian_f,trian_x)
    nothing
  elseif have_compatible_domains(trian_f,get_background_triangulation(trian_x))
    nothing
  elseif have_compatible_domains(get_background_triangulation(trian_f),get_background_triangulation(trian_x))
    nothing
  elseif have_compatible_domains(trian_x,get_background_triangulation(trian_f))
    @unreachable """\n
    CellField objects defined on a sub-triangulation cannot be integrated
    with a CellQuadrature defined on the underlying background mesh.

    This happens e.g. when trying to integrate a CellField defined on a Neumann boundary
    with a CellQuadrature defined on the underlying background mesh.
    """
  else
    @unreachable """\n
    Your are trying to integrate a CellField using a CellQuadrature defined on incompatible
    triangulations. Verify that either the two objects are defined in the same triangulation
    or that the triangulaiton of the CellField is the background triangulation of the CellQuadrature.
    """
  end

  b = change_domain(f,quad.trian,quad.domain_style)
  x = get_cell_points(quad)
  bx = b(x)
  if quad.domain_style == PhysicalDomain()
    lazy_map(IntegrationMap(),bx,quad.cell_weight)
  else
    cell_map = get_cell_map(quad.trian)
    cell_Jt = lazy_map(∇,cell_map)
    cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)
    lazy_map(IntegrationMap(),bx,quad.cell_weight,cell_Jtx)
  end
end

function integrate(a,quad::CellQuadrature)
  b = CellField(a,quad.trian,quad.domain_style)
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
"""
Contributions added to the cells of the background Triangulation.
"""
function get_cell_measure(trian::Triangulation)
  quad = CellQuadrature(trian,0)
  cell_to_dV = integrate(1,quad)
  cell_to_bgcell = get_cell_to_bgcell(trian)
  bgtrian = get_background_triangulation(trian)
  bgcell_to_dV = zeros(num_cells(bgtrian))
  _meas_K_fill!(bgcell_to_dV,cell_to_dV,cell_to_bgcell)
  bgcell_to_dV
end

function _meas_K_fill!(bgcell_to_dV,cell_to_dV,cell_to_bgcell)
  cache = array_cache(cell_to_dV)
  for (cell, bgcell) in enumerate(cell_to_bgcell)
    dV = getindex!(cache,cell_to_dV,cell)
    bgcell_to_dV[bgcell] += dV
  end
end
