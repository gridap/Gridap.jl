"""

The exported names are
$(EXPORTS)
"""
module CellData

using Test
using DocStringExtensions
using FillArrays

using Gridap.Helpers
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Integration

using Gridap.Arrays: get_array
using Gridap.Arrays: evaluate!

import Gridap.Integration: get_coordinates
import Gridap.Integration: get_weights

"""
Trait that signals if a CellData type is implemented in the physical or the reference domain
"""
abstract type DomainStyle end
struct ReferenceDomain <: DomainStyle end
struct PhysicalDomain <: DomainStyle end
DomainStyle(::T) where T = DomainStyle(T)

"""
Data associated with the cells of a Triangulation.
CellData objects behave as if they are defined in the physical space of the triangulation.
But in some cases they are implemented as reference quantities plus some transformation to the physical domain.
"""
abstract type CellData end

"""
Get the stored array of cell-wise data. It can be defined in the physical or the reference domain.
"""
get_cell_data(a::CellData) = @abstractmethod

"""
Tell if the stored array is in the reference or physical domain
"""
DomainStyle(::Type{<:CellData}) = @abstractmethod

"""
Return the underlying Triangulation object
"""
get_triangulation(a::CellData) = @abstractmethod

"""
Change the underlying data to the target domain
"""
change_domain(a::CellField,target_domain::DomainStyle) = change_domain(a,DomainStyle(a),target_domain)
change_domain(a::CellField,input_domain::T,target_domain::T) where T<: DomainStyle = a
change_domain(a::CellField,input_domain::DomainStyle,target_domain::DomainStyle) = @abstractmethod

# Tester
"""
"""
function test_cell_data(a::CellData)
  @test isa(get_cell_data(a),AbstractArray)
  @test isa(get_triangulation(a),Triangulation)
  @test isa(DomainStyle(a),DomainStyle)
end

# Some API

"""
Get the raw array of cell data defined in the physical space.
"""
get_array(a::CellData) = get_cell_data(change_domain(a,PhysicalDomain()))

"""
"""
Base.length(a::CellData) = length(get_cell_data(a))

# Evaluation / integration points

"""
A single point or an array of points on the cells of a Triangulation
"""
struct CellPoint{DS} <: CellData
  cell_point::AbstractArray{<:Union{Point,AbstractArray{<:Point}}}
  trian::Triangulation
  domain_style::DS
end

get_cell_data(f::CellPoint) = f.cell_point
get_triangulation(f::CellPoint) = f.trian
DomainStyle(::Type{CellPoint{DS}}) where DS = DS()

function change_domain(a::CellPoint,::ReferenceDomain,::PhysicalDomain)
  cell_map = get_cell_map(a.trian)
  cell_q = a.cell_point
  cell_x = lazy_map(evaluate,cell_map,cell_q)
  CellPoint(cell_x,a.trian,PhysicalDomain())
end

function change_domain(a::CellPoint,::PhysicalDomain,::ReferenceDomain)
  cell_map = get_cell_map(a.trian)
  cell_invmap = lazy_map(inverse_map,cell_map)
  cell_x = a.cell_point
  cell_q = lazy_map(evaluate,cell_invmap,cell_x)
  CellPoint(cell_q,a.trian,ReferenceDomain())
end

# Quadrature weights

"""
"""
struct CellWeight{DS} <: CellData
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
struct CellQuadrature{DS} <: CellData
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

get_cell_data(f::CellWeight) = f.cell_quad
get_triangulation(f::CellWeight) = f.trian
DomainStyle(::Type{CellWeight{DS}}) where DS = DS()

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

abstract type CellField <: CellData end

function change_domain(a::CellField,::ReferenceDomain,::PhysicalDomain)
  trian = get_triangulation(a)
  cell_map = get_cell_map(trian)
  cell_invmap = lazy_map(inverse_map,cell_map)
  cell_field_ref = get_cell_data(cell_field)
  cell_field_phys = lazy_map(Broadcasting(∘),cell_field_ref,cell_invmap)
  GenericCellField(cell_field_phys,trian,PhysicalDomain())
end

function change_domain(a::CellField,::PhysicalDomain,::ReferenceDomain)
  trian = get_triangulation(a)
  cell_map = get_cell_map(trian)
  cell_field_phys = get_cell_data(cell_field)
  cell_field_ref = lazy_map(Broadcasting(∘),cell_field_phys,cell_invmap)
  GenericCellField(cell_field_ref,trian,ReferenceDomain())
end

function change_domain(a::CellField,target_trian::Triangulation,target_domain::DomainStyle)
  change_domain(a,DomainStyle(a),target_trian,target_domain)
end

function change_domain(a::CellField,::ReferenceDomain,trian::Triangulation,::ReferenceDomain)
  trian_a = get_triangulation(a)
  if trian_a === trian
    return a
  elseif trian_a === get_background_triangulation(trian)
    cell_id = get_cell_id(trian)
    @notimplementedif isa(cell_id,SkeletonPair)
    cell_a_q = lazy_map(Reindex(get_cell_data(a)),cell_id)
    cell_s2q = get_cell_ref_map(trian)
    cell_field = lazy_map(Broadcasting(∘),cell_a_q,cell_s2q)
    GenericCellField(cell_field,trian,ReferenceDomain())
  else
    @unreachable """\n
    We cannot move the given CellField to the reference domain of the requested triangulation.
    Make sure that the given triangulation is either the same as the trianulation on which the
    CellField is defined, or that the latter triangulation is the background of the former.
    """
  end
end

function change_domain(a::CellField,::PhysicalDomain,trian::Triangulation,::PhysicalDomain)
  trian_a = get_triangulation(a)
  if trian_a === trian
    return a
  elseif trian_a === get_background_triangulation(trian)
    cell_id = get_cell_id(trian)
    @notimplementedif isa(cell_id,SkeletonPair)
    cell_field = lazy_map(Reindex(get_cell_data(a)),cell_id)
    GenericCellField(cell_field,trian,PhysicalDomain())
  else
    @unreachable """\n
    We cannot move the given CellField to the physical domain of the requested triangulation.
    Make sure that the given triangulation is either the same as the trianulation on which the
    CellField is defined, or that the latter triangulation is the background of the former.
    """
  end
end

function change_domain(a::CellField,::PhysicalDomain,trian::Triangulation,::ReferenceDomain)
  a_trian = change_domain(a,trian,PhysicalDomain())
  change_domain(a_trian,ReferenceDomain())
end

function change_domain(a::CellField,::ReferenceDomain,trian::Triangulation,::PhysicalDomain)
  a_phys = change_domain(a,PhysicalDomain())
  change_domain(a_phys,trian,PhysicalDomain())
end

struct GenericCellField{DS} <: CellField
  cell_field::AbstractArray{<:Union{Field,AbstractArray{<:Field}}}
  trian::Triangulation
  domain_style::DS
end

get_cell_data(f::GenericCellField) = f.cell_point
get_triangulation(f::GenericCellField) = f.trian
DomainStyle(::Type{GenericCellField{DS}}) where DS = DS()

# Evaluation of CellFields

function evaluate!(cache,f::CellField,x::Point)
  @notimplemented """\n
  Evaluation of CellField objects in a given Point is not yet implemented.

  This is a feature that we want to have at some point in Gridap.
  If you are ready to help with this implementation, please contact with the
  Gridap administrators.
  """
end

function evaluate!(cache,f::CellField,x::CellPoint)
  _f, _x = _to_common_domain(f,x)
  cell_field = get_cell_data(_f)
  cell_point = get_cell_data(_x)
  lazy_map(evaluate,cell_field,cell_point)
end

function _to_common_domain(f::CellField,x::CellPoint)

  trian_f = get_triangulation(f)
  trian_x = get_triangulation(x)

  if trian_f === trian_x
    nothing
  elseif trian_f === get_background_triangulation(trian_x)
    nothing
  elseif trian_x === get_background_triangulation(trian_f)
    @unreachable """\n
    CellField objects defined on a sub-triangulation cannot be evaluated
    on the underlying background mesh.
    """
  elseif
    @unreachable """\n
    Your are trying to evaluate a CellField on a CellPoint object defined on incompatible
    triangulations. Verify that either the two objects are defined in the same triangulation
    or that the triangulaiton of the CellField is the background triangulation of the CellPoint.
    """
  end

  f_on_trian_x = change_domain(f,trian_x,DomainStyle(x))
  f_on_trian_x, x
end

# Operations between CellField

#function evaluate!(cache,k::Operation,a::CellField,b::CellField)
#  a1, b1 = _to_common_domain(a,b)
#  @assert DomainStyle(a1) == DomainStyle(b1)
#  @assert get_triangulation(a1) == get_triangulation(b1)
#  trian = get_triangulation(a1)
#  domain_style = DomainStyle(a1)
#  cell_a = get_cell_data(a1)
#  cell_b = get_cell_data(a2)
#  cell_c = lazy_map(Broadcasting(k),cell_a,cell_b)
#  GenericCellField(cell_c,trian,domain_style)
#end

function evaluate!(cache,k::Operation,a::CellField...)

  b = _to_common_domain(a...)
  trian = get_triangulation(first(b))
  domain_style = DomainStyle(first(b))
  cell_b = map(get_cell_data,b)

  # Try to catch errors before hiding them deep in the lazy operation tree.
  if num_cells(trian) != 0
    bi = map(first,cell_b)
    ci = Broadcasting(k)(bi...)
    if domain_style == ReferenceDomain()
      xi = first(get_cell_ref_coordinates(trian))
    else
      xi = first(get_cell_coordinates(trian))
    end
    cixi = evaluate(ci,xi)
  end

  cell_c = lazy_map(Broadcasting(k),cell_b...)
  GenericCellField(cell_c,trian,domain_style)
end

#function _to_common_domain(a::CellField,b::CellField)
#  trian_a = get_triangulation(a)
#  trian_b = get_triangulation(b)
#  target_domain =  DomainStyle(a) == DomainStyle(b) ? DomainStyle(a) : ReferenceDomain()
#  if trian_a === trian_b
#    target_trian = trian_a
#  elseif trian_a === get_background_triangulation(trian_b)
#    target_trian = trian_b
#  elseif trian_b === get_background_triangulation(trian_a)
#    target_trian = trian_a
#  else
#    @unreachable """\n
#    You are trying to do an operation between CellField objects defined on
#    incompatible triangulations.
#
#    Make sure that both CellField objects are defined on the same triangulation
#    or that one triangulation is the background of the other.
#    """
#  end
#  a1 = change_domain(a,target_trian,target_domain)
#  b1 = change_domain(b,target_trian,target_domain)
#  a1, b1
#end

function _to_common_domain(a::CellField...)

  # Find a suitable domain style
  if any( map(i->DomainStyle(a)==ReferenceDomain(),a) )
    target_domain = ReferenceDomain()
  else
    target_domain = PhysicalDomain()
  end

  # Find a suitable triangulation
  msg = """\n
  You are trying to do an operation between CellField objects defined on
  incompatible triangulations.

  Make sure that all CellField objects are defined on the background triangulation
  or that the number of different sub-triangulations is equal to one.

  Examples:

  - Incompatible case: 2 cell fields defined on 2 different Neumann boundaries.

  - Compatible case: 3 cell fields 2, two them on the same Neumann boundary and the other on the background mesh.
  """
  trian_candidates = unique(objectid,map(get_triangulation,a))
  if length(trian_candidates) == 1
    target_trian = first(trian_candidates)
  elseif length(trian_candidates) == 2
    trian_a, trian_b = trian_candidates
    if trian_a === trian_b
      target_trian = trian_a
    elseif trian_a === get_background_triangulation(trian_b)
      target_trian = trian_b
    elseif trian_b === get_background_triangulation(trian_a)
      target_trian = trian_a
    else
      @unreachable msg
    end
  else
    @unreachable msg
  end
  map(i->change_domain(i,target_trian,target_domain),a)
end

end # module
