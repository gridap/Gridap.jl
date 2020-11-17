
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

"""
"""
function CellQuadrature(trian::Triangulation,degree::Integer)
  ctype_to_reffe, cell_to_ctype = compress_cell_data(get_cell_reffe(trian))
  ctype_to_quad = map(r->Quadrature(get_polytope(r),degree),ctype_to_reffe)
  ctype_to_point = map(get_coordinates,ctype_to_quad)
  ctype_to_weigth = map(get_weights,ctype_to_quad)
  cell_quad = expand_cell_data(ctype_to_quad,cell_to_ctype)
  cell_point = expand_cell_data(ctype_to_point,cell_to_ctype)
  cell_weight = expand_cell_data(ctype_to_weigth,cell_to_ctype)
  CellQuadrature(cell_quad,cell_point,cell_weight,trian,ReferenceDomain())
end

function CellQuadrature(trian::AppendedTriangulation,degree1::Integer,degree2::Integer)
  quad1 = CellQuadrature(trian.a,degree1)
  quad2 = CellQuadrature(trian.b,degree2)
  lazy_append(quad1,quad2,trian)
end

function CellQuadrature(trian::AppendedTriangulation,degree::Integer)
  CellQuadrature(trian,degree,degree)
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

get_cell_data(f::CellQuadrature) = f.cell_quad
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



