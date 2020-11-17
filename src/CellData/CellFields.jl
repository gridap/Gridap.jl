
"""
A single point or an array of points on the cells of a Triangulation
CellField objects can be evaluated efficiently at CellPoint instances.
"""
struct CellPoint{DS} <: CellDatum
  cell_ref_point::AbstractArray{<:Union{Point,AbstractArray{<:Point}}}
  cell_phys_point::AbstractArray{<:Union{Point,AbstractArray{<:Point}}}
  trian::Triangulation
  domain_style::DS
end

function CellPoint(
  cell_ref_point::AbstractArray{<:Union{Point,AbstractArray{<:Point}}},
  trian::Triangulation,
  domain_style::ReferenceDomain)

  cell_map = get_cell_map(trian)
  cell_phys_point = lazy_map(evaluate,cell_map,cell_ref_point)
  CellPoint(cell_ref_point,cell_phys_point,trian,domain_style)
end

function CellPoint(
  cell_phys_point::AbstractArray{<:Union{Point,AbstractArray{<:Point}}},
  trian::Triangulation,
  domain_style::PhysicalDomain)
  cell_invmap = lazy_map(inverse_map,cell_map)

  cell_map = get_cell_map(trian)
  cell_invmap = lazy_map(inverse_map,cell_map)
  cell_ref_point = lazy_map(evaluate,cell_invmap,cell_phys_point)
  CellPoint(cell_ref_point,cell_phys_point,trian,domain_style)
end

function get_cell_data(f::CellPoint)
  if DomainStyle(f) == ReferenceDomain()
    f.cell_ref_point
  else
    f.cell_phys_point
  end
end

get_triangulation(f::CellPoint) = f.trian
DomainStyle(::Type{CellPoint{DS}}) where DS = DS()

function change_domain(a::CellPoint,::ReferenceDomain,::PhysicalDomain)
  CellPoint(a.cell_ref_point,a.cell_phys_point,a.trian,PhysicalDomain())
end

function change_domain(a::CellPoint,::PhysicalDomain,::ReferenceDomain)
  CellPoint(a.cell_ref_point,a.cell_phys_point,a.trian,ReferenceDomain())
end

# Possibly with a different name
"""
"""
function get_cell_points(trian::Triangulation)
  cell_ref_coords = get_cell_ref_coordinates(trian)
  cell_phys_coords = get_cell_coordinates(trian)
  CellPoint(cell_ref_coords,cell_phys_coords,trian,ReferenceDomain())
end

"""
"""
abstract type CellField <: CellDatum end

function CellField(f::Function,trian::Triangulation,domain_style::DomainStyle)
  s = size(get_cell_map(trian))
  cell_field = Fill(GenericField(f),s)
  GenericCellField(cell_field,trian,PhysicalDomain())
end

function CellField(f::Number,trian::Triangulation,domain_style::DomainStyle)
  s = size(get_cell_map(trian))
  cell_field = Fill(ConstantField(f),s)
  GenericCellField(cell_field,trian,domain_style)
end

function CellField(f::AbstractArray{<:Number},trian::Triangulation,domain_style::DomainStyle)
  @check length(f)==num_cells(trian)  """\n
  You are trying to build a CellField from an array of length $(length(f))
  on a Triangulation with $(num_cells(trian)) cells. The length of the given array
  and the number of cells should match.
  """
  cell_field = lazy_map(ConstantField,f)
  GenericCellField(cell_field,trian,domain_style)
end

function CellField(f::CellField,trian::Triangulation,domain_style::DomainStyle)
  change_domain(f,trian,domain_style)
end

function CellField(f,trian::Triangulation)
  CellField(f,trian,ReferenceDomain())
end

function get_normal_vector(trian::Triangulation)
  cell_normal = get_facet_normal(trian)
  @assert ! isa(cell_normal,SkeletonPair)
  GenericCellField(cell_normal,trian,ReferenceDomain())
end

evaluate!(cache,f::Function,x::CellPoint) = CellField(f,get_triangulation(x))(x)

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
  cell_field_phys = get_cell_data(a)
  cell_field_ref = lazy_map(Broadcasting(∘),cell_field_phys,cell_map)
  GenericCellField(cell_field_ref,trian,ReferenceDomain())
end

"""
"""
function change_domain(a::CellField,target_trian::Triangulation,target_domain::DomainStyle)
  change_domain(a,DomainStyle(a),target_trian,target_domain)
end

function change_domain(a::CellField,::ReferenceDomain,trian::Triangulation,::ReferenceDomain)
  trian_a = get_triangulation(a)
  if have_compatible_domains(trian_a,trian)
    return a
  elseif have_compatible_domains(trian_a,get_background_triangulation(trian))
    cell_id = get_cell_id(trian)
    @assert ! isa(cell_id,SkeletonPair)
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
  if have_compatible_domains(trian_a,trian)
    return a
  elseif have_compatible_domains(trian_a,get_background_triangulation(trian))
    cell_id = get_cell_id(trian)
    @assert ! isa(cell_id,SkeletonPair)
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

"""
"""
struct GenericCellField{DS} <: CellField
  cell_field::AbstractArray{<:Union{Field,AbstractArray{<:Field}}}
  trian::Triangulation
  domain_style::DS
end

get_cell_data(f::GenericCellField) = f.cell_field
get_triangulation(f::GenericCellField) = f.trian
DomainStyle(::Type{GenericCellField{DS}}) where DS = DS()

# Evaluation of CellFields

(a::CellField)(x) = evaluate(a,x)

function evaluate!(cache,f::CellField,x::Point)
  @notimplemented """\n
  Evaluation of a CellField at a given Point is not implemented yet.

  This is a feature that we want to have at some point in Gridap.
  If you are ready to help with this implementation, please contact the
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

  if have_compatible_domains(trian_f,trian_x)
    nothing
  elseif have_compatible_domains(trian_f,get_background_triangulation(trian_x))
    nothing
  elseif have_compatible_domains(trian_x,get_background_triangulation(trian_f))
    @unreachable """\n
    CellField objects defined on a sub-triangulation cannot be evaluated
    on the underlying background mesh.

    This happens e.g. when trying to evaluate a CellField defined on a Neumann boundary
    at a CellPoint defined on the underlying background mesh.
    """
  else
    @unreachable """\n
    Your are trying to evaluate a CellField on a CellPoint object defined on incompatible
    triangulations. Verify that either the two objects are defined in the same triangulation
    or that the triangulaiton of the CellField is the background triangulation of the CellPoint.
    """
  end

  f_on_trian_x = change_domain(f,trian_x,DomainStyle(x))
  f_on_trian_x, x
end

# Gradient

function gradient(a::CellField)
  cell_∇a = lazy_map(Broadcasting(∇),get_cell_data(a))
  if DomainStyle(a) == PhysicalDomain()
    g = cell_∇a
  else
    cell_map = get_cell_map(get_triangulation(a))
    g = lazy_map(Broadcasting(push_∇),cell_∇a,cell_map)
  end
  GenericCellField(g,get_triangulation(a),DomainStyle(a))
end

function ∇∇(a::CellField)
  cell_∇∇a = lazy_map(Broadcasting(∇∇),get_cell_data(a))
  if DomainStyle(a) == PhysicalDomain()
    h = cell_∇∇a
  else
    cell_map = get_cell_map(get_triangulation(a))
    h = lazy_map(Broadcasting(push_∇∇),cell_∇∇a,cell_map)
  end
  GenericCellField(h,get_triangulation(a),DomainStyle(a))
end

# Operations between CellField

function evaluate!(cache,k::Operation,a::CellField...)
  _operate_cellfields(k,a...)
end

function evaluate!(cache,k::Operation,a::Union{Function,CellField}...)
  _operate_cellfields(k,_convert_to_cellfields(a...)...)
end

function evaluate!(cache,k::Operation,a::Union{Number,CellField}...)
  _operate_cellfields(k,_convert_to_cellfields(a...)...)
end

function evaluate!(cache,k::Operation,a::Union{AbstractArray{<:Number},CellField}...)
  _operate_cellfields(k,_convert_to_cellfields(a...)...)
end

# Why julia hangs with this method????
#
#function evaluate!(cache,k::Operation,a::Union{Function,Number,AbstractArray,CellField}...)
#  b = _convert_to_cellfields(a...)
#  _operate_cellfields(k,b...)
#end

struct OperationCellField{DS} <: CellField
  op::Operation
  args::Tuple
  trian::Triangulation
  domain_style::DS
  function OperationCellField(op::Operation,args::CellField...)

    @assert length(args) > 0
    trian = get_triangulation(first(args))
    domain_style = DomainStyle(first(args))
    @check all( map(i->DomainStyle(i)==domain_style,args) )
    @check all( map(i->have_compatible_domains(get_triangulation(i),trian),args) )

    if num_cells(trian)>0
      x = get_cell_points(trian)
      try
         ax = map(i->i(x),args)
         axi = map(first,ax)
         r = BroadcastingFieldOpMap(op.op)(axi...)
      catch
        @unreachable """\n
        It is not possible to perform operation $(op.op) on the given cell fields.

        See the catched error for more information.
        """
      end
    end

    new{typeof(domain_style)}(op,args,trian,domain_style)
  end
end

function get_cell_data(f::OperationCellField)
  a = map(get_cell_data,f.args)
  lazy_map(Broadcasting(f.op),a...)
end
get_triangulation(f::OperationCellField) = f.trian
DomainStyle(::Type{OperationCellField{DS}}) where DS = DS()

function evaluate!(cache,f::OperationCellField,x::CellPoint)
  ax = map(i->i(x),f.args)
  lazy_map(BroadcastingFieldOpMap(f.op.op),ax...)
end

function change_domain(f::OperationCellField,target_trian::Triangulation,target_domain::DomainStyle)
  args = map(i->change_domain(i,target_trian,target_domain),f.args)
  OperationCellField(f.op,args...)
end

function _operate_cellfields(k::Operation,a...)
  b = _to_common_domain(a...)
  OperationCellField(k,b...)
end

function _convert_to_cellfields(a...)
  a1 = filter(i->isa(i,CellField),a)
  a2 = _to_common_domain(a1...)
  target_domain = DomainStyle(first(a2))
  target_trian = get_triangulation(first(a2))
  map(i->CellField(i,target_trian,target_domain),a)
end

function _to_common_domain(a::CellField...)

  # Find a suitable domain style
  if any( map(i->DomainStyle(i)==ReferenceDomain(),a) )
    target_domain = ReferenceDomain()
  else
    target_domain = PhysicalDomain()
  end

  # Find a suitable triangulation
  msg = """\n
  You are trying to operate CellField objects defined on incompatible triangulations.

  Make sure that all CellField objects are defined on the background triangulation
  or that the number of different sub-triangulations is equal to one.

  For instace:

  - 3 cell fields 2, two them on the same Neumann boundary and the other on the background mesh is OK.

  - 2 cell fields defined on 2 different Neumann boundaries is NOT OK.
  """
  trian_candidates = unique(objectid,map(get_triangulation,a))
  if length(trian_candidates) == 1
    target_trian = first(trian_candidates)
  elseif length(trian_candidates) == 2
    trian_a, trian_b = trian_candidates
    if have_compatible_domains(trian_a,trian_b)
      target_trian = trian_a
    elseif have_compatible_domains(trian_a,get_background_triangulation(trian_b))
      target_trian = trian_b
    elseif have_compatible_domains(trian_b,get_background_triangulation(trian_a))
      target_trian = trian_a
    else
      @unreachable msg
    end
  else
    @unreachable msg
  end
  map(i->change_domain(i,target_trian,target_domain),a)
end

# Composition (this replaces the @law macro)

Base.:(∘)(f::Function,g::CellField) = Operation(f)(g)
Base.:(∘)(f::Function,g::Tuple{Vararg{CellField}}) = Operation(f)(g...)
Base.:(∘)(f::Function,g::Tuple{Vararg{Union{AbstractArray{<:Number},CellField}}}) = Operation(f)(g...)
Base.:(∘)(f::Function,g::Tuple{Vararg{Union{Function,CellField}}}) = Operation(f)(g...)

# Define some of the well known arithmetic ops

# Unary ops

for op in (:symmetric_part,:inv,:det,:abs,:abs2,:+,:-,:tr,:transpose,:adjoint,:grad2curl)
  @eval begin
    ($op)(a::CellField) = Operation($op)(a)
  end
end

# Binary ops

for op in (:inner,:outer,:double_contraction,:+,:-,:*,:cross,:dot,:/)
  @eval begin
    ($op)(a::CellField,b::CellField) = Operation($op)(a,b)
    ($op)(a::CellField,b::Number) = Operation($op)(a,b)
    ($op)(a::Number,b::CellField) = Operation($op)(a,b)
    ($op)(a::CellField,b::Function) = Operation($op)(a,b)
    ($op)(a::Function,b::CellField) = Operation($op)(a,b)
    ($op)(a::CellField,b::AbstractArray{<:Number}) = Operation($op)(a,b)
    ($op)(a::AbstractArray{<:Number},b::CellField) = Operation($op)(a,b)
  end
end

# Skeleton related Operations

function Base.getproperty(x::CellField, sym::Symbol)
  if sym in (:⁺,:plus)
    CellFieldAt{:plus}(x)
  elseif sym in (:⁻, :minus)
    CellFieldAt{:minus}(x)
  else
    getfield(x, sym)
  end
end

function Base.propertynames(x::CellField, private=false)
  (fieldnames(typeof(x))...,:⁺,:plus,:⁻,:minus)
end

struct CellFieldAt{T,F} <: CellField
  parent::F
  CellFieldAt{T}(parent::CellField) where T = new{T,typeof(parent)}(parent)
end

get_cell_data(f::CellFieldAt) = get_cell_data(f.parent)
get_triangulation(f::CellFieldAt) = get_triangulation(f.parent)
DomainStyle(::Type{CellFieldAt{T,F}}) where {T,F} = DomainStyle(F)
gradient(a::CellFieldAt{P}) where P = CellFieldAt{P}(gradient(a.parent))
∇∇(a::CellFieldAt{P}) where P = CellFieldAt{P}(∇∇(a.parent))

function CellFieldAt{T}(parent::OperationCellField) where T
  args = map(i->CellFieldAt{T}(i),parent.args)
  OperationCellField(parent.op,args...)
end

function get_normal_vector(trian::SkeletonTriangulation)
  cell_normal_plus = get_facet_normal(trian.plus)
  cell_normal_minus = get_facet_normal(trian.minus)
  plus = GenericCellField(cell_normal_plus,trian,ReferenceDomain())
  minus = GenericCellField(cell_normal_minus,trian,ReferenceDomain())
  SkeletonPair(plus,minus)
end

for op in (:outer,:*,:dot)
  @eval begin
    ($op)(a::CellField,b::SkeletonPair{<:CellField}) = Operation($op)(a,b)
    ($op)(a::SkeletonPair{<:CellField},b::CellField) = Operation($op)(a,b)
  end
end

function evaluate!(cache,k::Operation,a::CellField,b::SkeletonPair{<:CellField})
  plus = k(a.plus,b.plus)
  minus = k(a.minus,b.minus)
  SkeletonPair(plus,minus)
end

function evaluate!(cache,k::Operation,a::SkeletonPair{<:CellField},b::CellField)
  plus = k(a.plus,b.plus)
  minus = k(a.minus,b.minus)
  SkeletonPair(plus,minus)
end

jump(a::CellField) = a.⁺ - a.⁻
jump(a::SkeletonPair{<:CellField}) = a.⁺ + a.⁻ # a.⁻ results from multiplying by n.⁻. Thus we need to sum.

mean(a::CellField) = Operation(_mean)(a.⁺,a.⁻)
_mean(x,y) = 0.5*x + 0.5*y

# This is the fundamental part to make operations on the skeleton work.

function change_domain(a::CellField,target_trian::SkeletonTriangulation,target_domain::DomainStyle)
  trian_a = get_triangulation(a)
  if have_compatible_domains(trian_a,target_trian)
    return a
  elseif have_compatible_domains(trian_a,get_background_triangulation(target_trian))
    # In this case, we can safely take either plus or minus arbitrarily.
    if isa(a,GenericCellField) && (isa(a.cell_field,Fill{<:ConstantField}) || isa(a.cell_field,Fill{<:GenericField{<:Function}}))
      a_on_target_trian = change_domain(a,target_trian.plus,target_domain)
      return GenericCellField(get_cell_data(a_on_target_trian),target_trian,target_domain)
    else
      @unreachable """\n
      It is not possible to use the given CellField on a SkeletonTriangulation.
      Make sure that you are specifying which of the two possible traces,
      either plus (aka ⁺) or minus (aka ⁻) you want to use.
      """
    end
  else
    @unreachable """\n
    We cannot move the given CellField to the requested triangulation.
    Make sure that the given CellField is defined on the triangulation you want to work with.
    """
  end
end

function change_domain(a::CellFieldAt,trian::SkeletonTriangulation,target_domain::DomainStyle)
  trian_a = get_triangulation(a)
  if have_compatible_domains(trian_a,get_background_triangulation(trian))
    plus, minus = change_domain_skeleton(a.parent,trian,target_domain)
    if isa(a,CellFieldAt{:plus})
      return plus
    elseif isa(a,CellFieldAt{:minus})
      return minus
    else
      @unreachable
    end
  else
    @unreachable """\n
    It is not allowd to writte `u.⁺` of `u.⁻` for the given CellField.
    Make sure that the CellField `u` is either defined on the background mesh
    or it is a normal vector extracted from a SkeletonTriangulation.
    """
  end
end

function change_domain_skeleton(a::CellField,trian::SkeletonTriangulation,target_domain::DomainStyle)
  a_on_plus_trian = change_domain(a,trian.plus,target_domain)
  a_on_minus_trian = change_domain(a,trian.minus,target_domain)
  plus = GenericCellField(get_cell_data(a_on_plus_trian),trian,target_domain)
  minus = GenericCellField(get_cell_data(a_on_minus_trian),trian,target_domain)
  plus, minus
end

function change_domain(f::OperationCellField,target_trian::SkeletonTriangulation,target_domain::DomainStyle)
  args = map(i->change_domain(i,target_trian,target_domain),f.args)
  OperationCellField(f.op,args...)
end

# Just to provide more meaningful error messages
function (a::SkeletonPair{<:CellField})(x)
  @unreachable """\n
  You are trying to evaluate a CellField on a mesh skeleton but you have not specified which of the
  two sides i.e. plus (aka ⁺) or minus (aka ⁻) you want to select.

  For instance, if you have extracted the normal vector and the cell points from a SkeletonTriangulation

      x = get_cell_points(strian)
      n = get_normal_vector(strian)

  Evaluating `n(x)` is not allowed. You need to call either `n.⁺(x)` or `n.⁻(x)`.
  """
end
