"""
A single point or an array of points belonging to cells of a Triangulation.
`CellField`s can be evaluated efficiently at `CellPoint`s.
"""
struct CellPoint{DS,A,B,C} <: CellDatum
  cell_ref_point::A
  cell_phys_point::B
  trian::C
  domain_style::DS
end

function CellPoint(
  cell_ref_point::AbstractArray,
  trian::Triangulation,
  domain_style::ReferenceDomain)

  cell_map = get_cell_map(trian)
  cell_phys_point = lazy_map(evaluate,cell_map,cell_ref_point)
  CellPoint(cell_ref_point,cell_phys_point,trian,domain_style)
end

function CellPoint(
  cell_phys_point::AbstractArray,
  trian::Triangulation,
  domain_style::PhysicalDomain)
  cell_map = get_cell_map(trian)
  cell_invmap = lazy_map(inverse_map,cell_map)
  cell_ref_point = lazy_map(evaluate,cell_invmap,cell_phys_point)
  CellPoint(cell_ref_point,cell_phys_point,trian,domain_style)
end

function get_data(f::CellPoint)
  if DomainStyle(f) == ReferenceDomain()
    f.cell_ref_point
  else
    f.cell_phys_point
  end
end

get_triangulation(f::CellPoint) = f.trian
DomainStyle(::Type{CellPoint{DS,A,B,C}}) where {DS,A,B,C} = DS()

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

function Base.:(==)(a::CellPoint,b::CellPoint)
  a.trian == b.trian &&
  a.cell_ref_point == b.cell_ref_point &&
  a.cell_phys_point == b.cell_phys_point &&
  a.domain_style == b.domain_style
end

"""
    abstract type CellField <: CellDatum end

Represents a `Field` that is defined cell-wise, that is in each cell, of a
`Triangulation`.

Can only be efficiently evaluated at [`CellPoint`](@ref),
since otherwise a cell map inversion is required to find the cell owning the
point one wishes to evaluate the cell field at.
"""
abstract type CellField <: CellDatum end

function similar_cell_field(f::CellField,cell_data,trian,ds)
  GenericCellField(cell_data,trian,ds)
end

@inline function similar_cell_field(f,cell_data)
  similar_cell_field(f,cell_data,get_triangulation(f),DomainStyle(f))
end

function Base.show(io::IO,::MIME"text/plain",f::CellField)
  show(io,f)
  print(io,":")
  print(io,"\n num_cells: $(num_cells(f))")
  print(io,"\n DomainStyle: $(DomainStyle(f))")
  print(io,"\n Triangulation: $(get_triangulation(f))")
  print(io,"\n Triangulation id: $(objectid(get_triangulation(f)))")
end


"""
    CellField(f::Function, trian::Triangulation, ds::DomainStyle=PhysicalDomain())
    CellField(f::Number,   trian::Triangulation, ds::DomainStyle)
    CellField(f::AbstractArray{<:Number}, trian::Triangulation, ds::DomainStyle)
    CellField(f::CellField, trian::Triangulation, ds::DomainStyle=ReferenceDomain())

Generic `CellField` factory for closed form cell fields. The default `DomainStyle` is `ReferenceDomain`
except for julia `Function`s.
""" # Functions are assumed physical
CellField(f::Function, trian::Triangulation) = CellField(f, trian, PhysicalDomain())

function CellField(f::Function, trian::Triangulation, ::PhysicalDomain)
  s = size(get_cell_map(trian))
  cell_field = Fill(GenericField(f), s)
  GenericCellField(cell_field, trian, PhysicalDomain())
end

function CellField(f::Function, trian::Triangulation, ::ReferenceDomain)
  phys_f = CellField(f, trian, PhysicalDomain())
  change_domain(phys_f, PhysicalDomain(), ReferenceDomain())
end

function CellField(f::Number,trian::Triangulation,domain_style::DomainStyle)
  s = size(get_cell_map(trian))
  cell_field = Fill(constant_field(f),s)
  GenericCellField(cell_field,trian,domain_style)
end

function CellField(f::AbstractArray{<:Number},trian::Triangulation,domain_style::DomainStyle)
  @check length(f)==num_cells(trian)  """\n
  You are trying to build a CellField from an array of length $(length(f))
  on a Triangulation with $(num_cells(trian)) cells. The length of the given array
  and the number of cells should match.
  """
  cell_field = lazy_map(constant_field,f)
  GenericCellField(cell_field,trian,domain_style)
end

function CellField(f::CellField,trian::Triangulation,domain_style::DomainStyle)
  change_domain(f,trian,domain_style)
end

function CellField(f,trian::Triangulation)
  CellField(f,trian,ReferenceDomain())
end

"""
    get_normal_vector(trian::Triangulation)
"""
function get_normal_vector(trian::Triangulation)
  cell_normal = get_facet_normal(trian)
  get_normal_vector(trian, cell_normal)
end

"""
    get_tangent_vector(trian::Triangulation)
"""
function get_tangent_vector(trian::Triangulation)
  cell_tangent = get_edge_tangent(trian)
  get_tangent_vector(trian, cell_tangent)
end

function get_normal_vector(trian::Triangulation,cell_vectors::AbstractArray)
  GenericCellField(cell_vectors,trian,ReferenceDomain())
end

function get_tangent_vector(trian::Triangulation,cell_vectors::AbstractArray)
  GenericCellField(cell_vectors,trian,ReferenceDomain())
end

evaluate!(cache,f::Function,x::CellPoint) = CellField(f,get_triangulation(x))(x)

function change_domain(a::CellField,::ReferenceDomain,::PhysicalDomain)
  trian = get_triangulation(a)
  cell_map = get_cell_map(trian)
  cell_invmap = lazy_map(inverse_map,cell_map)
  cell_field_ref = get_data(a)
  cell_field_phys = lazy_map(Broadcasting(∘),cell_field_ref,cell_invmap)
  similar_cell_field(a,cell_field_phys,trian,PhysicalDomain())
end

function change_domain(a::CellField,::PhysicalDomain,::ReferenceDomain)
  trian = get_triangulation(a)
  cell_map = get_cell_map(trian)
  cell_field_phys = get_data(a)
  cell_field_ref = lazy_map(Broadcasting(∘),cell_field_phys,cell_map)
  similar_cell_field(a,cell_field_ref,trian,ReferenceDomain())
end

function change_domain(a::CellField,target_trian::Triangulation,target_domain::DomainStyle)
  change_domain(a,get_triangulation(a),DomainStyle(a),target_trian,target_domain)
end

function change_domain(a::CellField,source_domain::DomainStyle,target_trian::Triangulation,target_domain::DomainStyle)
  change_domain(a,get_triangulation(a),source_domain,target_trian,target_domain)
end

function change_domain(a::CellField,strian::Triangulation,::ReferenceDomain,ttrian::Triangulation,::ReferenceDomain)
  msg = """\n
  We cannot move the given CellField to the reference domain of the requested triangulation.
  Make sure that the given triangulation is either the same as the triangulation on which the
  CellField is defined, or that the latter triangulation is the background of the former.
  """
  if strian === ttrian
    return a
  end
  @check is_change_possible(strian,ttrian) msg
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  change_domain_ref_ref(a,ttrian,sglue,tglue)
end

function change_domain(a::CellField,strian::Triangulation,::PhysicalDomain,ttrian::Triangulation,::PhysicalDomain)
  msg = """\n
  We cannot move the given CellField to the physical domain of the requested triangulation.
  Make sure that the given triangulation is either the same as the triangulation on which the
  CellField is defined, or that the latter triangulation is the background of the former.
  """
  if strian === ttrian
    return a
  end
  @check is_change_possible(strian,ttrian) msg
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  change_domain_phys_phys(a,ttrian,sglue,tglue)
end

function change_domain(a::CellField,strian::Triangulation,::PhysicalDomain,trian::Triangulation,::ReferenceDomain)
  a_trian = change_domain(a,trian,PhysicalDomain())
  change_domain(a_trian,ReferenceDomain())
end

function change_domain(a::CellField,strian::Triangulation,::ReferenceDomain,trian::Triangulation,::PhysicalDomain)
  a_phys = change_domain(a,PhysicalDomain())
  change_domain(a_phys,trian,PhysicalDomain())
end

function change_domain_ref_ref(
  a::CellField,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::FaceToFaceGlue)
  sface_to_field = get_data(a)
  mface_to_sface = sglue.mface_to_tface
  tface_to_mface = tglue.tface_to_mface
  tface_to_mface_map = tglue.tface_to_mface_map
  mface_to_field = extend(sface_to_field,mface_to_sface)
  tface_to_field_s = lazy_map(Reindex(mface_to_field),tface_to_mface)
  tface_to_field_t = lazy_map(Broadcasting(∘),tface_to_field_s,tface_to_mface_map)
  similar_cell_field(a,tface_to_field_t,ttrian,ReferenceDomain())
end

function change_domain_phys_phys(
  a::CellField,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::FaceToFaceGlue)
  sface_to_field = get_data(a)
  mface_to_sface = sglue.mface_to_tface
  tface_to_mface = tglue.tface_to_mface
  mface_to_field = extend(sface_to_field,mface_to_sface)
  tface_to_field = lazy_map(Reindex(mface_to_field),tface_to_mface)
  similar_cell_field(a,tface_to_field,ttrian,PhysicalDomain())
end

"""
"""
struct GenericCellField{DS} <: CellField
  cell_field::AbstractArray
  trian::Triangulation
  domain_style::DS
  function GenericCellField(
    cell_field::AbstractArray,
    trian::Triangulation,
    domain_style::DomainStyle)

    DS = typeof(domain_style)
    new{DS}(Fields.MemoArray(cell_field),trian,domain_style)
  end
end

get_data(f::GenericCellField) = f.cell_field
get_triangulation(f::GenericCellField) = f.trian
DomainStyle(::Type{GenericCellField{DS}}) where DS = DS()

(a::CellField)(x) = evaluate(a,x)

function evaluate!(cache,f::CellField,x::CellPoint)
  _f, _x = _to_common_domain(f,x)
  cell_field = get_data(_f)
  cell_point = get_data(_x)
  lazy_map(evaluate,cell_field,cell_point)
end

function _to_common_domain(f::CellField,x::CellPoint)
  trian_f = get_triangulation(f)
  trian_x = get_triangulation(x)
  f_on_trian_x = change_domain(f,trian_x,DomainStyle(x))
  f_on_trian_x, x
end

# Differentiation operators
#
# ∂ᴺ of a CellField `a` returns the field whose values are the N-th physical
# derivative, expressed in `a`'s own coordinates, its DomainStyle. The low
# level differentiation implementation differentiates in `a`'s coordinates and
# converts the derivation indices from `a`'s coordinate basis to the physical
# one.

function gradient(a::CellField)
  cell_∇a = lazy_map(Broadcasting(∇),get_data(a))
  if DomainStyle(a) == PhysicalDomain()
    g = cell_∇a
  else
    cell_map = get_cell_map(get_triangulation(a))
    g = lazy_map(Broadcasting(push_∇),cell_∇a,cell_map)
  end
  similar_cell_field(a,g,get_triangulation(a),DomainStyle(a))
end

function (k::PartialDerivatives{N})(a::CellField) where N
  cell_∂a = lazy_map(Broadcasting(k), get_data(a))
  if DomainStyle(a) == PhysicalDomain()
    ka = cell_∂a
  else
    @notimplementedif N != 1 """\n
    Partial derivatives of order $N of a cell field defined in the reference
    domain are not implemented yet: their pushforward by the cell map involves
    the derivatives of the map itself, up to order $N, as `push_∇∇` does.
    """
    cell_map = get_cell_map(get_triangulation(a))
    cell_Jt = lazy_map(partial_derivatives, cell_map)
    cell_iJt = lazy_map(Operation(pinvJt), cell_Jt)
    # `pushforward_partials` transports the derivation index with `pinv(Jt)`
    ka = lazy_map(Broadcasting(Operation(push_first_partials)), cell_∂a, cell_iJt)
  end
  similar_cell_field(a, ka, get_triangulation(a), DomainStyle(a))
end

function DIV(a::CellField)
  DIVa = DIV(get_data(a))
  if DomainStyle(a) == PhysicalDomain()
    @notimplemented
  end
  similar_cell_field(a,DIVa,get_triangulation(a),DomainStyle(a))
end

function ∇∇(a::CellField)
  cell_∇∇a = lazy_map(Broadcasting(∇∇),get_data(a))
  if DomainStyle(a) == PhysicalDomain()
    h = cell_∇∇a
  else
    cell_map = get_cell_map(get_triangulation(a))
    # this assumes that cell_map is affine, and push_∇∇ enforces it.
    h = lazy_map(Broadcasting(push_∇∇),cell_∇∇a,cell_map)
  end
  similar_cell_field(a,h,get_triangulation(a),DomainStyle(a))
end

# This function has to be removed when ∇⋅∇(a) is implemented
laplacian(a::CellField) = tr(∇∇(a))

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

struct OperationCellField{DS,O} <: CellField
  op::Operation{O}
  args::Tuple
  trian::Triangulation
  domain_style::DS
  memo::Dict{Any,Any}
  function OperationCellField(op::Operation,args::CellField...)

    @assert length(args) > 0
    trian = get_triangulation(first(args))
    domain_style = DomainStyle(first(args))
    @check all(i -> DomainStyle(i) == domain_style, args)

    # This is only to catch errors in user code
    # as soon as possible.
    # if num_cells(trian) > 0
    #   @check begin
    #     pts = _get_cell_points(args...)
    #     #x = testitem(get_data(pts))
    #     #f = map(ak -> testitem(get_data(ak)), args)
    #     #fx = map(fk -> return_value(fk,x), f)
    #     #r = Fields.BroadcastingFieldOpMap(op.op)(fx...)
    #     ax = map(i->i(pts),args)
    #     axi = map(first,ax)
    #     r = Fields.BroadcastingFieldOpMap(op.op)(axi...)
    #     true
    #   end
    # end

    DS, O = typeof(domain_style), typeof(op.op)
    new{DS,O}(op,args,trian,domain_style,Dict())
  end
end

function _get_cell_points(args::CellField...)
  k = findfirst(i->isa(i,CellState),args)
  if k === nothing
    j = findall(i->isa(i,OperationCellField),args)
    if length(j) == 0
      _get_cell_points(first(args))
    else
      _get_cell_points(args[j]...)
    end
  else
    args[k].points
  end
end

function _get_cell_points(a::CellField)
  trian = get_triangulation(a)
  get_cell_points(trian)
end

function _get_cell_points(a::OperationCellField...)
  b = []
  for ai in a
    for i in ai.args
      push!(b,i)
    end
  end
  _get_cell_points(b...)
end

function _get_cell_points(a::OperationCellField)
  _get_cell_points(a.args...)
end

function get_data(f::OperationCellField)
  a = map(get_data,f.args)
  lazy_map(Broadcasting(f.op),a...)
end
get_triangulation(f::OperationCellField) = f.trian
DomainStyle(::Type{OperationCellField{DS,O}}) where {DS,O} = DS()

function evaluate!(cache,f::OperationCellField,x::CellPoint)
  ax = map(i->i(x),f.args)
  lazy_map(Fields.BroadcastingFieldOpMap(f.op.op),ax...)
end

for diffop in (:gradient,:divergence,:curl,:∇∇)
  for op in (:+,:-)
    # Sometimes we cannot call `get_data` from the `OperatorCellField`. We dispatch
    # sooner so that we can evaluate the gradient in specific cases.
    # This is necessary to take the gradient of the sum of two FEBasis, for instance.
    @eval begin
      function Fields.$diffop(f::OperationCellField{DS,typeof($op)}) where {DS}
        args = map($diffop,f.args)
        OperationCellField(Operation($op),args...)
      end
    end
  end
end

for op in (:+,:-)
  @eval begin
    function (k::PartialDerivatives)(f::OperationCellField{DS,typeof($op)}) where DS
      args = map(k, f.args)
      OperationCellField(Operation($op), args...)
    end
  end
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

  For instance:

  - 3 cell fields 2, two them on the same Neumann boundary and the other on the background mesh is OK.

  - 2 cell fields defined on 2 different Neumann boundaries is NOT OK.
  """
  trian_candidates = unique(objectid,map(get_triangulation,a))
  if length(trian_candidates) == 1
    target_trian = first(trian_candidates)
  elseif length(trian_candidates) == 2
    trian_a, trian_b = trian_candidates
    sa_tb = is_change_possible(trian_a,trian_b)
    sb_ta = is_change_possible(trian_b,trian_a)
    if sa_tb && sb_ta
      target_trian = best_target(trian_a,trian_b)
    elseif !sa_tb && sb_ta
      target_trian = trian_a
    elseif sa_tb && !sb_ta
      target_trian = trian_b
    else
      @unreachable msg
    end
  else
    m = """\n
    Cannote operate cellfields defined over more than 2 different
    triangulations at this moment.
    """
    @notimplemented
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

for op in (:symmetric_part,:inv,:det,:abs,:abs2,:+,:-,:tr,:transpose,:adjoint,
  :grad2curl,:real,:imag,:conj,:flat_hodge_star,:flat_flat,:flat_sharp,
  :flat_vol_coeff,:from_0form,:to_1form,:from_1form
)
  @eval begin
    ($op)(a::CellField) = Operation($op)(a)
  end
end

# Binary ops

for op in (:outer,:double_contraction,:+,:-,:*,:cross,:dot,:/,
  :interior_product, :∧, :flat_inner, :flat_alt_inner, :pullback, :pushforward
)
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

# Metric dependent ops: the triangulation's metric is automatically fed to the
# operator, but an explicitly given is not changed.

function _operate_with_metric(op, euclidean_op, a...)
  a_cf = _convert_to_cellfields(a...)
  # The metric is that of the triangulation the operation is evaluated on, so
  # it is read after the arguments are brought to their common domain.
  md_cf = get_metric_data(get_triangulation(first(a_cf)))
  _metric_operation(op, euclidean_op, md_cf, a_cf...)
end

_metric_operation(op, euclidean_op, ::Tuple{UniformScaling}, a::CellField...) =
  Operation(euclidean_op)(a...)
_metric_operation(op, euclidean_op, metric_data::Tuple, a::CellField...) =
  Operation(op)(a..., _metric_args(op, metric_data)...)

# metric_data is `(g, g_inv, sqrt_abs_g)`, this eagerly computes all metric args
# the operation might require
_metric_args(::typeof(hodge_star), metric_data::Tuple) = metric_data
_metric_args(::typeof(flat),       metric_data::Tuple) = (metric_data[1],)
_metric_args(::typeof(sharp),      metric_data::Tuple) = (metric_data[2],)
_metric_args(::typeof(inner),      metric_data::Tuple) = metric_data[1:2]
_metric_args(::typeof(alt_inner),  metric_data::Tuple) = metric_data[1:2]

# One field argument
for (op, euclidean_op) in (
  (:hodge_star,:flat_hodge_star), (:flat,:flat_flat), (:sharp,:flat_sharp)
)
  @eval begin
    ($op)(a::CellField) = _operate_with_metric($op,$euclidean_op,a)
    ($op)(a::CellField, ::UniformScaling) = Operation($euclidean_op)(a)
    ($op)(a::CellField, g, gs...) = Operation($op)(a,g,gs...)
  end
end

# Two field arguments
for (op, euclidean_op) in ((:inner,:flat_inner), (:alt_inner,:flat_alt_inner))
  for (Ta,Tb) in (
    (:CellField,:CellField),
    (:CellField,:Number), (:Number,:CellField),
    (:CellField,:Function), (:Function,:CellField),
    (:CellField,:(AbstractArray{<:Number})), (:(AbstractArray{<:Number}),:CellField)
  )
    @eval ($op)(a::$Ta,b::$Tb) = _operate_with_metric($op,$euclidean_op,a,b)
  end
  @eval begin
    ($op)(a::CellField, b::CellField, ::UniformScaling) = Operation($euclidean_op)(a,b)
    ($op)(a::CellField, b::CellField, g, gs...) = Operation($op)(a, b, g, gs...)
  end
end

Base.broadcasted(f,a::CellField,b::CellField) = Operation((i,j)->f.(i,j))(a,b)
Base.broadcasted(f,a::Number,b::CellField) = Operation((i,j)->f.(i,j))(a,b)
Base.broadcasted(f,a::CellField,b::Number) = Operation((i,j)->f.(i,j))(a,b)
Base.broadcasted(f,a::Function,b::CellField) = Operation((i,j)->f.(i,j))(a,b)
Base.broadcasted(f,a::CellField,b::Function) = Operation((i,j)->f.(i,j))(a,b)
Base.broadcasted(::typeof(*),::typeof(∇),f::CellField) = Operation(Fields._extract_grad_diag)(∇(f))
Base.broadcasted(::typeof(*),s::Fields.ShiftedNabla,f::CellField) = Operation(Fields._extract_grad_diag)(s(f))

dot(::typeof(∇),f::CellField) = divergence(f)
function (*)(::typeof(∇),f::CellField)
  msg = "Syntax ∇*f has been removed, use ∇⋅f (\\nabla \\cdot f) instead"
  error(msg)
end
outer(::typeof(∇),f::CellField) = gradient(f)
outer(f::CellField,::typeof(∇)) = transpose(gradient(f))
cross(::typeof(∇),f::CellField) = curl(f)

"""
    get_physical_coordinate(trian::Triangulation)

In contrast to get_cell_map, the returned object:
- is a [`CellField`](@ref)
- its gradient is the identity tensor
"""
function get_physical_coordinate(trian::Triangulation)
  CellField(_phys_coord,trian)
end

_phys_coord(x) = x

_phys_coord_grad(x) = one(typeof(outer(x,x)))

gradient(::typeof(_phys_coord)) = _phys_coord_grad

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

function Base.propertynames(x::CellField, private::Bool=false)
  (fieldnames(typeof(x))...,:⁺,:plus,:⁻,:minus)
end

struct CellFieldAt{T,F} <: CellField
  parent::F
  CellFieldAt{T}(parent::CellField) where T = new{T,typeof(parent)}(parent)
end

get_data(f::CellFieldAt) = get_data(f.parent)
get_triangulation(f::CellFieldAt) = get_triangulation(f.parent)
DomainStyle(::Type{CellFieldAt{T,F}}) where {T,F} = DomainStyle(F)
gradient(a::CellFieldAt{P}) where P = CellFieldAt{P}(gradient(a.parent))
∇∇(a::CellFieldAt{P}) where P = CellFieldAt{P}(∇∇(a.parent))
(k::PartialDerivatives)(a::CellFieldAt{P}) where P = CellFieldAt{P}(k(a.parent))
function similar_cell_field(f::CellFieldAt{T},cell_data,trian,ds) where T
  parent = similar_cell_field(f.parent,cell_data,trian,ds)
  CellFieldAt{T}(parent)
end

function CellFieldAt{T}(parent::OperationCellField) where T
  args = map(i->CellFieldAt{T}(i),parent.args)
  OperationCellField(parent.op,args...)
end

function get_normal_vector(trian::Triangulation,cell_normal::SkeletonPair)
  plus = get_normal_vector(trian,cell_normal.plus)
  minus = get_normal_vector(trian,cell_normal.minus)
  SkeletonPair(plus,minus)
end

function get_tangent_vector(trian::Triangulation,cell_tangent::SkeletonPair)
  plus = get_tangent_vector(trian,cell_tangent.plus)
  minus = get_tangent_vector(trian,cell_tangent.minus)
  SkeletonPair(plus,minus)
end

for op in (:outer,:*,:dot,:inner,:/)
  @eval begin
    ($op)(a::CellField,b::SkeletonPair{<:CellField}) = Operation($op)(a,b)
    ($op)(a::SkeletonPair{<:CellField},b::CellField) = Operation($op)(a,b)
    ($op)(a::SkeletonPair{<:CellField},b::SkeletonPair{<:CellField}) = Operation($op)(a,b)
  end
end

function evaluate!(cache,k::Operation,a::SkeletonPair{<:CellField})
  plus = k(a.plus)
  minus = k(a.minus)
  SkeletonPair(plus,minus)
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

function evaluate!(cache,k::Operation,a::SkeletonPair{<:CellField},b::SkeletonPair{<:CellField})
  plus = k(a.plus,b.plus)
  minus = k(a.minus,b.minus)
  SkeletonPair(plus,minus)
end

"""
    jump(a::CellField)
    jump(a::SkeletonPair{<:CellField})

Jump operator at interior facets of the supporting `Triangulation`, defined by
`jump`(`a` n) = ⟦`a` n⟧ = `a`⁺n⁺ + `a`⁻n⁻, where n is an oriented normal field
to the interior facets, n⁺ = -n⁻ are the normal pointing into the element on the +
and - side of the facets, and `a`⁺/`a`⁻ are the restrictions of `a` to each
element respectively.
"""
jump(a::CellField) = a.⁺ - a.⁻
jump(a::SkeletonPair{<:CellField}) = a.⁺ + a.⁻ # a.⁻ results from multiplying by n.⁻. Thus we need to sum.

"""
    mean(a::CellField)

Similar to [`jump`](@ref), but for the mean operator `a` ⟶ (`a`⁺ + `a`⁻)/2.
"""
mean(a::CellField) = Operation(_mean)(a.⁺,a.⁻)
_mean(x,y) = 0.5*x + 0.5*y

# This is the fundamental part to make operations on the skeleton work.

for fun in (:change_domain_ref_ref,:change_domain_phys_phys)
  @eval begin

    function $fun(
      a::CellField,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::SkeletonPair)
      msg = """\n
        It is not possible to use the given CellField on a SkeletonTriangulation.
        Make sure that you are specifying which of the two possible traces,
        either plus (aka ⁺) or minus (aka ⁻) you want to use.
        """
      # If the underlying array is Fill it does not matter if we restrict if from
      # the plus or minus side.
      @check isa(get_array(a.cell_field),Fill)
      plus = $fun(a,ttrian,sglue,tglue.plus)
      plus
    end

    function $fun(
      a::CellFieldAt,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::SkeletonPair)
      if isa(a,CellFieldAt{:plus})
        $fun(a,ttrian,sglue,tglue.plus)
      elseif isa(a,CellFieldAt{:minus})
        $fun(a,ttrian,sglue,tglue.minus)
      else
        @unreachable
      end
    end

    function $fun(
      f::OperationCellField,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::SkeletonPair)
      args = map(i->$fun(i,ttrian,sglue,tglue),f.args)
      OperationCellField(f.op,args...)
    end

  end
end

function change_domain(a::SkeletonPair, ::ReferenceDomain, ::PhysicalDomain)
  plus = change_domain(a.plus,ReferenceDomain(),PhysicalDomain())
  minus = change_domain(a.minus,ReferenceDomain(),PhysicalDomain())
  return SkeletonPair(plus,minus)
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

################################################################################
# Differential geometry operators needing derivatives or evaluation point

function exterior_derivative(a::CellField)
  cell_a = get_data(a)
  if DomainStyle(a) == ReferenceDomain() && is_pushed_by(cell_a, TensorPushforward())
    # `𝑑` does commute with the pullback, so on a basis pushed forward by
    # `TensorPushforward` it is taken in the reference space and pushed forward
    # as a `k+1` form.
    cell_da = exterior_derivative(cell_a)
    # where exterior_derivative is overloaded for field arrays, in ReferenceFEs/Pullbacks.jl
  else
    # A field that is not a pushforward for reference space (e.g. physical or
    # 0-form reference) has it's partial derivatives computed in the physical
    # space.
    ∂a = partial_derivatives(a)
    cell_da = lazy_map(Broadcasting(Operation(partials_to_exteriorderivative)), get_data(∂a))
  end
  similar_cell_field(a, cell_da, get_triangulation(a), DomainStyle(a))
end

# `δ` does not comute with pullback, the derivatives have to be taken in physical space
function codifferential(a::CellField)
  trian = get_triangulation(a)
  codifferential(a, get_metric_data(trian), get_metric_partials_data(trian))
end

codifferential(a::CellField, ::UniformScaling) =
  codifferential(a, (LinearAlgebra.I,), LinearAlgebra.I)

function codifferential(a::CellField, ::Tuple{UniformScaling}, ::UniformScaling)
  ∂a = partial_derivatives(a)
  cell_δa = lazy_map(Broadcasting(Operation(partials_to_flat_codifferential)), get_data(∂a))
  similar_cell_field(a, cell_δa, get_triangulation(a), DomainStyle(a))
end

function codifferential(a::CellField, metric_data::Tuple, ∂g::CellField)
  g_inv = metric_data[2]
  Γ = Operation(partials_to_Christoffel_2)(∂g, g_inv)
  Operation(partials_to_codifferential)(partial_derivatives(a), a, Γ, g_inv)
end

# `κ_x(a) = ι_x(a)` contracts `a` with the position vector
koszul(a::CellField) = Operation(koszul)(get_physical_coordinate(get_triangulation(a)), a)

