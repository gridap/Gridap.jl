
"""
"""
struct CellWeight{DS} <: CellDatum
  cell_point::AbstractArray{<:AbstractArray{<:Point}}
  cell_weight::AbstractArray{<:AbstractArray{<:Real}}
  trian::Triangulation
  domain_style::DS
end

get_cell_data(f::CellWeight) = f.cell_weight
get_triangulation(f::CellWeight) = f.trian
DomainStyle(::Type{CellWeight{DS}}) where DS = DS()

function change_domain(a::CellWeight,::ReferenceDomain,::PhysicalDomain)
  cell_map = get_cell_map(a.trian)
  cell_Jt = lazy_map(∇,cell_map)
  cell_det = lazy_map(Operation(det),cell_Jt)
  cell_q = a.cell_point
  cell_det_q = lazy_map(evaluate,cell_det,cell_q)
  cell_weight = lazy_map(Broadcasting(*),cell_det_q,a.cell_weight)
  cell_x = lazy_map(evaluate,cell_map,cell_q)
  CellWeight(cell_x,cell_weight,a.trian,PhysicalDomain())
end

function change_domain(a::CellWeight,::PhysicalDomain,::ReferenceDomain)
  @notimplemented
end

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
  ctype_to_reffe = get_reffes(trian)
  cell_to_ctype = get_cell_type(trian)
  ctype_to_quad = map(r->Quadrature(get_polytope(r),degree),ctype_to_reffe)
  ctype_to_point = map(get_coordinates,ctype_to_quad)
  ctype_to_weigth = map(get_coordinates,ctype_to_quad)
  cell_point = expand_cell_data(ctype_to_point,cell_to_ctype)
  cell_weight = expand_cell_data(ctype_to_weigth,cell_to_ctype)
  CellQuadrature(cell_quad,cell_point,cell_weight,trian,ReferenceDomain())
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

function get_coordinates(a::CellQuadrature)
  CellPoint(a.cell_point,a.trian,a.domain_style)
end

function get_weights(a::CellQuadrature)
  CellWeight(a.cell_weight,a.trian,a.domain_style)
end

# 
