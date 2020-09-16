
# DISCLAIMER
# ==========
#  This is a experimental implementation of the new behavior of CellFields
#  Most of this functionality can (should) be improved and implemented in more low level
#  modules. The issue is mainly code style not performance.

"""
    abstract type CellField <: GridapType end

This is a plain array of fields plus some metadata.
It is interpreted as a `Field` object on each cell of a computational
mesh.
"""
abstract type CellField <: GridapType end

"""
    get_array(cf::CellField)
"""
function get_array(cf::CellField)
  @abstractmethod
end

"""
    get_cell_axes(cf::CellField)
"""
function get_cell_axes(cf::CellField)
  @abstractmethod
end

"""
    get_memo(cf::CellField)
"""
function get_memo(cf::CellField)
  @abstractmethod
end

"""
This trait provides information about the field size at a single
evaluation point.

For physical fields `MetaSizeStyle(T) == Val( () )`
For test bases `MetaSizeStyle(T) == Val( (:,) )`
For trial bases `MetaSizeStyle(T) = Val( (1,:) )`
For fields representing elemental matrices `MetaSizeStyle(T) = Val( (:,:) )`
"""
MetaSizeStyle(::Type{<:CellField}) = @abstractmethod
MetaSizeStyle(::T) where T = MetaSizeStyle(T)
get_metasize(::Type{T}) where T = get_val_parameter(MetaSizeStyle(T))
get_metasize(::T) where T = get_val_parameter(MetaSizeStyle(T))
is_test(::Type{T}) where T = get_metasize(T) == (:,)
is_test(::T) where T = is_test(T)
is_trial(::Type{T}) where T = get_metasize(T) == (1,:)
is_trial(::T) where T = is_trial(T)
is_basis(::Type{T}) where T =  is_test(T) || is_trial(T)
is_basis(::T) where T = is_basis(T)

function Base.getproperty(x::CellField, sym::Symbol)
  if sym in (:inward,:⁺)
    get_inward(x)
  elseif sym in (:outward,:⁻)
    get_outward(x)
  else
    getfield(x, sym)
  end
end

# Tester

"""
    test_cell_field(
      cf::CellField,
      x::AbstractArray,
      b::AbstractArray,
      pred=(==);
      grad=nothing)
"""
function test_cell_field(cf::CellField,x::CellPoint,b::AbstractArray,pred=(==);grad=nothing)
  a = evaluate(cf,x)
  test_array(a,b,pred)
  if grad != nothing
    g = evaluate(gradient(cf),x)
    test_array(g,grad,pred)
  end
  @test get_metasize(cf) == get_metasize(cf)
  @test get_cell_axes(cf) == get_cell_axes(cf)
  @test isa(get_cell_axes(cf),AbstractArray{<:Tuple})
  @test isa(get_metasize(cf),Tuple)
  @test isa(MetaSizeStyle(cf),Val)
  @test isa(get_memo(cf),Dict)
end

function test_cell_field(cf::CellField,x::AbstractArray,b::AbstractArray,pred=(==);grad=nothing)
  test_cell_field(cf,GenericCellPoint(x),b,pred;grad=grad)
end

# Concrete implementation

"""
struct GenericCellField{S} <: CellField
  # Private fields
end

"""
struct GenericCellField{S} <: CellField
  array::AbstractArray
  cell_axes::AbstractArray
  metasize::Val{S}
  memo::Dict
  function GenericCellField(
    array::AbstractArray,
    cell_axes::AbstractArray,
    metasize::Val{S}) where S
    memo = Dict()
    new{S}(array,cell_axes,metasize,memo)
  end
end

function GenericCellField(array::AbstractArray)
  cell_axes = Fill((),length(array))
  metasize = Val(())
  GenericCellField(array,cell_axes,metasize)
end

function get_array(cf::GenericCellField)
  cf.array
end

function get_cell_axes(cf::GenericCellField)
  cf.cell_axes
end

function get_memo(cf::GenericCellField)
  cf.memo
end

function MetaSizeStyle(::Type{<:GenericCellField{S}}) where {S}
  Val{S}()
end

# The rest of the file is some default API using the CellField interface

function get_wrapped_cell_field(cf::CellField)
  cf
end

function wrap_cell_field(cf::CellField,a::CellField)
  a
end

function align_cell_fields(cfs::CellField...)
  cf1 = first(cfs)
  map(cf -> align_cell_fields(cf1,cf)[2], cfs)
end

function align_cell_fields(a::CellField)
  a
end

function align_cell_fields(a::CellField,b::CellField)
  map(c->GenericCellField(get_array(c),get_cell_axes(c),MetaSizeStyle(c)), (a,b))
end

function align_cell_fields(a::GenericCellField,b::GenericCellField)
  a,b
end

"""
    Base.:(==)(a::CellField,b::CellField)

Check if two `CellField` objects are the same. Uses `===` by default. 
"""
function Base.:(==)(a::CellField,b::CellField)
  if length(a) != length(b)
    return false
  end
  a === b || get_array(a) === get_array(b) || get_array(a) == get_array(b)
end

"""
    length(cf::CellField)
"""
function Base.length(cf::CellField)
  length(get_cell_axes(cf))
end

"""
    evaluate(f::CellField,x::CellPoint)

It is assumed that `x` is in the domain of the array of fields
in `get_array(f)`.
"""
function evaluate(f::CellField,x)
  key = (:evaluate,objectid(x))
  memo = get_memo(f)
  if !haskey(memo,key)
    memo[key] = compute_evaluate(f,x)
  end
  memo[key]
end

function compute_evaluate(f::CellField,x::CellPoint)
  b = reindex(get_array(f),get_cell_id(x))
  evaluate_field_array(b,get_array(x))
end

function compute_evaluate(f::CellField,x::AbstractArray)
  compute_evaluate(f,GenericCellPoint(x))
end

"""
    (f::CellField)(x)

Functor like-evaluation `f(x)` for `CellField` objects. Syntactic sugar for `evaluate(f,x)`.
"""
function (f::CellField)(x)
  evaluate(f,x)
end

"""
    Base.:∘(f::CellField,ϕ::CellMap)

Composition `(f∘ϕ)(q) == f(ϕ(q)) `
"""
function Base.:∘(f::CellField,ϕ::CellMap)
  _compose_cell_fields(f,ϕ)
end

function _compose_cell_fields(f,ϕ)
  g = reindex(f,get_cell_id(ϕ))
  array = compose_field_arrays(get_array(g),get_array(ϕ))
  GenericCellField(array,get_cell_axes(g),MetaSizeStyle(g))
end

function Base.:∘(object,ϕ::CellMap)
  f = convert_to_cell_field(object,num_cell_ids(ϕ))
  f∘ϕ
end

Base.:∘(f::CellField,ϕ::AppendedCellMap) = lazy_append(f∘ϕ.ϕ1,f∘ϕ.ϕ2)

Base.:∘(f,ϕ::AppendedCellMap) = lazy_append(f∘ϕ.ϕ1,f∘ϕ.ϕ2)

function Arrays.lazy_append(a::CellField,b::CellField)
  array = lazy_append(get_array(a),get_array(b))
  cell_axs = lazy_append(get_cell_axes(a),get_cell_axes(b))
  GenericCellField(array,cell_axs,MetaSizeStyle(a))
end

"""
    reindex(f::CellField,a::AbstractArray)
"""
function Arrays.reindex(f::CellField,a::AbstractArray)
  array = reindex(get_array(f),a)
  cell_axes = reindex(get_cell_axes(f),a)
  GenericCellField(array,cell_axes,MetaSizeStyle(f))
end

# Bases-related
function trialize_cell_basis(test::CellField)
  _trialize_cell_basis(test,MetaSizeStyle(test))
end

function _trialize_cell_basis(test,metasize)
  @unreachable
end

function _trialize_cell_basis(trial,metasize::Val{(1,:)})
  trial
end

function _trialize_cell_basis(test::CellField,metasize::Val{(:,)})
  a = get_wrapped_cell_field(test)
  array = trialize_array_of_bases(get_array(a))
  axs = apply(Fields._add_singleton_block,get_cell_axes(a))
  b = GenericCellField(array,axs,Val((1,:)))
  wrap_cell_field(test,b)
end

function Fields.lincomb(cf::CellField,b::AbstractArray)
  a = get_wrapped_cell_field(cf)
  @notimplementedif !is_test(a)
  @assert length(a) == length(b)
  array = lincomb(get_array(a),b)
  axs = Fill((),length(b))
  c = GenericCellField(array,axs,Val(()))
  wrap_cell_field(cf,c)
end

# Diff ops

"""
    gradient(cf::CellField)
"""
function gradient(cf::CellField)
  key = :gradient
  memo = get_memo(cf)
  if !haskey(memo,key)
    memo[key] = compute_gradient(cf)
  end
  memo[key]
end

function compute_gradient(cf::CellField)
  a = get_array(cf)
  ag = field_array_gradient(a)
  GenericCellField(ag,get_cell_axes(cf),MetaSizeStyle(cf))
end

# Operations

function operate(op,cf::CellField)
  key = objectid(op)
  memo = get_memo(cf)
  if !haskey(memo,key)
    a = get_wrapped_cell_field(cf)
    b = operate_arrays_of_fields(Fields._UnimplementedField,op,get_array(a))
    c = GenericCellField(b,get_cell_axes(a),MetaSizeStyle(a))
    memo[key] = wrap_cell_field(cf,c)
  end
  memo[key]
end

function operate(op,_cf1::CellField,_cf2::CellField)
  cf1, cf2 = align_cell_fields(_cf1,_cf2)
  @assert length(cf1) == length(cf2)
  a1 = get_wrapped_cell_field(cf1)
  a2 = get_wrapped_cell_field(cf2)
  b = operate_arrays_of_fields(Fields._UnimplementedField,op,get_array(a1),get_array(a2))
  axs = apply(field_operation_axes,get_cell_axes(a1),get_cell_axes(a2))
  metasize = field_operation_metasize(get_metasize(a1),get_metasize(a2))
  c = GenericCellField(b,axs,Val(metasize))
  wrap_cell_field(cf1,c)
end

function operate(op,cf1::CellField,object)
  cf2 = convert_to_cell_field(object,length(cf1))
  operate(op,cf1,cf2)
end

function operate(op,object,cf2::CellField)
  cf1 = convert_to_cell_field(object,length(cf2))
  operate(op,cf1,cf2)
end

function operate(op,_args::CellField...)
  _args_aligned = align_cell_fields(_args...)
  args = map(get_wrapped_cell_field, _args_aligned)
  arrs = map(get_array,args)
  m = operate_arrays_of_fields(Fields._UnimplementedField,op,arrs...)
  axs = apply(field_operation_axes,map(get_cell_axes,args)...)
  metasize = field_operation_metasize(map(get_metasize,args)...)
  c = GenericCellField(m,axs,Val(metasize))
  wrap_cell_field(first(_args_aligned),c)
end

function operate(op,a1::CellField,args...)
  cfs = map(obj->convert_to_cell_field(obj,length(a1)),args)
  operate(op,a1,cfs...)
end

# Conversions

function convert_to_cell_field(f,l::Integer)
  @abstractmethod
end

function convert_to_cell_field(f::CellField,l::Integer)
  @assert length(f) == l
  f
end

function convert_to_cell_field(f::AbstractArray{<:Number},l::Integer)
  @assert length(f) == l
  GenericCellField(f)
end

function convert_to_cell_field(f::Number,l::Integer)
  GenericCellField(Fill(f,l))
end

function convert_to_cell_field(f::Function,l::Integer)
  ff = function_field(f)
  GenericCellField(Fill(ff,l))
end

# Optimizations associated with inverse maps (i.e., with the reference space)

function Base.:∘(f::CellField,ϕinv::InverseCellMap)
  CellFieldComposedWithInverseMap(f,ϕinv)
end

function Base.:∘(f::CellMap,ϕinv::InverseCellMap)
  g = GenericCellField(get_array(f))
  g∘ϕinv
end

# TODO gradient not yet properly computed.
# We assume that the gradients of the given field are already in the 
# physical space.
"""
    struct CellFieldComposedWithInverseMap{T<:CellField} <: CellField
      f::T
      ϕinv::InverseCellMap
    end
struct representing `f∘inverse_map(ϕ)`, i.e. a `CellField` in the reference space
pushed to the physical one.  This struct allows a number of optimizations.
In particular, the computation of `(f∘inverse_map(ϕ))∘ϕ`, which is simply computed as`f`.
It also computes `∇(f∘inverse_map(ϕ))` as `(∇(ϕ)\\∇(f))∘inverse_map(ϕ)` which is in turn
stored in an instance of `CellFieldComposedWithInverseMap`.
"""
struct CellFieldComposedWithInverseMap{T<:CellField} <: CellField
  f::T
  ϕinv::InverseCellMap
  memo::Dict
  function CellFieldComposedWithInverseMap(f::CellField,ϕinv::InverseCellMap)
    memo = Dict()
    T = typeof(f)
    new{T}(f,ϕinv,memo)
  end
end

Base.length(a::CellFieldComposedWithInverseMap) = length(a.ϕinv)

function get_array(a::CellFieldComposedWithInverseMap)
  get_array(a.f∘get_inverse_cell_map(a))
end

function get_wrapped_cell_field(a::CellFieldComposedWithInverseMap)
  a.f
end

function wrap_cell_field(a::CellFieldComposedWithInverseMap,b::CellField)
  CellFieldComposedWithInverseMap(b,a.ϕinv)
end

function get_cell_axes(a::CellFieldComposedWithInverseMap)
  @notimplementedif length(a.f) != length(a.ϕinv.ϕ)
  get_cell_axes(a.f)
end

MetaSizeStyle(::Type{CellFieldComposedWithInverseMap{T}}) where T = MetaSizeStyle(T)

get_memo(a::CellFieldComposedWithInverseMap) = a.memo

function get_inverse_cell_map(a::CellFieldComposedWithInverseMap)
  ϕinv = GenericCellMap(get_array(a.ϕinv),get_cell_id(a.ϕinv),num_cell_ids(a.ϕinv))
  ϕinv
end

function Arrays.reindex(f::CellFieldComposedWithInverseMap,a::AbstractVector)
  CellFieldComposedWithInverseMap(reindex(f.f,a),reindex(f.ϕinv,a))
end

function compute_gradient(f::CellFieldComposedWithInverseMap)
  ag = field_array_gradient(get_array(f.f))
  fg = GenericCellField(ag,get_cell_axes(f.f),MetaSizeStyle(f.f))
  wrap_cell_field(f,fg)
end

# TODO
##function _compute_phys_gradient(f,f_ref)
##  ∇ϕ = gradient(f.ϕinv.ϕ)
##  ∇f = gradient(f_ref)
##  g = _phys_grad(∇f,∇ϕ)
##  CellFieldComposedWithInverseMap(g,f.ϕinv)
##end
##
##function _compute_phys_gradient(f,cf::CellFieldFromOperation{typeof(lincomb)})
##  a = _compute_phys_gradient(f,cf.args[1])
##  b = CellFieldFromOperation(cf.op,lincomb,(a.f,),get_cell_axes(cf),MetaSizeStyle(cf))
##  c = GenericCellField(get_array(b))
##  CellFieldComposedWithInverseMap(c,f.ϕinv)
##end
##
##function _phys_grad(∇f,∇ϕ)
##  #inv(∇ϕ)⋅∇f This also works but larger op tree
##  a = apply_to_field_array(Fields.PhysGrad(),get_array(∇f),get_array(∇ϕ))
##  GenericCellField(a,get_cell_axes(∇f),MetaSizeStyle(∇f))
##end

function Base.:∘(f::CellFieldComposedWithInverseMap,ϕ::CellMap)
  if f.ϕinv.ϕ == ϕ
    return f.f
  else
    return (f.f∘get_inverse_cell_map(f))∘ϕ
  end
end

function maps_can_simplify(f,ϕ::CellMap)
  if f.ϕinv.ϕ == ϕ
    true
  else
    false
  end
end

function Base.:∘(f::CellFieldComposedWithInverseMap,ϕ::ReindexedCellMap)
  if f.ϕinv.ϕ == ϕ
    return f.f
  else
    return reindex(f∘ϕ.cell_map,ϕ.ids)
  end
end

function maps_can_simplify(f,ϕ::ReindexedCellMap)
  if f.ϕinv.ϕ == ϕ
    return true
  else
    return maps_can_simplify(f,ϕ.cell_map)
  end
end

function Base.:∘(f::CellFieldComposedWithInverseMap,ϕ::FaceMap)
  if f.ϕinv.ϕ == ϕ
    return f.f
  else
    return (f∘ϕ.cell_map)∘ϕ.refface_to_refcell_map
  end
end

function maps_can_simplify(f,ϕ::FaceMap)
  if f.ϕinv.ϕ == ϕ
    return true
  else
    return maps_can_simplify(f,ϕ.cell_map)
  end
end

function Base.:∘(f::CellFieldComposedWithInverseMap,ϕ::AppendedCellMap)
  if f.ϕinv.ϕ == ϕ
    return f.f
  else
    return lazy_append(f∘ϕ.ϕ1,f∘ϕ.ϕ2)
  end
end

function maps_can_simplify(f,ϕ::AppendedCellMap)
  if f.ϕinv.ϕ == ϕ
    return true
  else
    return maps_can_simplify(f,ϕ.ϕ1) && maps_can_simplify(f,ϕ.ϕ2)
  end
end

function align_cell_fields(a::CellFieldComposedWithInverseMap,b::CellField)
  _align_cell_fields_a_wins(a,b)
end

function align_cell_fields(a::CellField,b::CellFieldComposedWithInverseMap)
  _align_cell_fields_b_wins(a,b)
end

function align_cell_fields(a::CellFieldComposedWithInverseMap,b::CellFieldComposedWithInverseMap)
  if a.ϕinv.ϕ == b.ϕinv.ϕ
    return a,b
  elseif maps_can_simplify(a,b.ϕinv.ϕ)
    return _align_cell_fields_b_wins(a,b)
  elseif maps_can_simplify(b,a.ϕinv.ϕ)
    return _align_cell_fields_a_wins(a,b)
  else
    @notimplemented
  end
end

function _align_cell_fields_a_wins(a,b)
  c = b∘a.ϕinv.ϕ
  a, wrap_cell_field(a,c)
end

function _align_cell_fields_b_wins(a,b)
  c = a∘b.ϕinv.ϕ
  wrap_cell_field(b,c), b
end

# Skeleton related

struct SkeletonFaceMap <:CellMap
  left::CellMap
  right::CellMap
  memo::Dict
  function SkeletonFaceMap(left::CellMap,right::CellMap)
    new(left,right,Dict())
  end
  function SkeletonFaceMap(left::FaceMap,right::FaceMap)
    _right = change_face_map(right,left.face_map)
    new(left,_right,Dict())
  end
end

Arrays.get_array(a::SkeletonFaceMap) = get_array(a.left)
get_cell_id(a::SkeletonFaceMap) = get_cell_id(a.left)
num_cell_ids(a::SkeletonFaceMap) = num_cell_ids(a.left)
get_memo(a::SkeletonFaceMap) = a.memo

#struct SkeletonCellField{F<:CellField} <: CellField
#  f::F
#  memo::Dict
#  SkeletonCellField(f::CellField) = new{typeof(f)}(f,Dict())
#end
#
#Arrays.get_array(a::SkeletonCellField) = get_array(a.f)
#get_cell_axes(a::SkeletonCellField) = get_cell_axes(a.f)
#get_memo(a::SkeletonCellField) = a.memo
#MetaSizeStyle(::Type{<:SkeletonCellField{F}}) where F = MetaSizeStyle(F)

function Base.:∘(f::CellField,ϕinv::InverseCellMap{<:SkeletonFaceMap})
  a = CellFieldComposedWithInverseMap(f,ϕinv)
  CellFieldFromOperation(identity,(a,))
end

struct CellFieldFromOperation{F} <: CellField
  op
  args::Tuple
  memo::Dict
  function CellFieldFromOperation(op,::F,args::Tuple) where F
    new{F}(op,args,Dict())
  end
end

function CellFieldFromOperation(op,args::Tuple)
  CellFieldFromOperation(op,op,args)
end

get_memo(a::CellFieldFromOperation) = a.memo

function operate(op,cfs::CellFieldFromOperation...)
  _operate_cfs(op,cfs...)
end

function operate(op,a::CellFieldFromOperation,b::CellField)
  _operate_cfs(op,a,b)
end

function operate(op,a::CellField,b::CellFieldFromOperation)
  _operate_cfs(op,a,b)
end

function operate(op,a::CellFieldFromOperation,b)
  _operate_cfs(op,a,b)
end

function operate(op,a,b::CellFieldFromOperation)
  _operate_cfs(op,a,b)
end

function _operate_cfs(op,cfs...)
  f(args...) = operate(op,args...)
  CellFieldFromOperation(f,op,cfs)
end

for op in (:get_inward, :get_outward)
  @eval begin

    function ($op)(a::CellField)
      CellFieldFromOperation($op,(a,))
    end

    function ($op)(a::CellFieldFromOperation{typeof(identity)})
      a
    end

    function ($op)(f::CellFieldFromOperation)
      f.op(map($op,f.args)...)
    end

    function gradient(a::CellFieldFromOperation{typeof($op)})
      key = :gradient
      memo = get_memo(a)
      if !haskey(memo,key)
        memo[key] = $op(gradient(first(a.args)))
      end
      memo[key]
    end

  end
end

function gradient(a::CellFieldFromOperation)
   key = :gradient
   memo = get_memo(a)
   if !haskey(memo,key)
     memo[key] = CellFieldFromOperation(gradient,(a,))
   end
   memo[key]
end

function Base.:∘(f::CellFieldFromOperation,ϕ::CellMap)
  f.op(map(i->i∘ϕ,f.args)...)
end

function Base.:∘(f::CellFieldFromOperation{typeof(get_inward)},ϕ::CellMap)
  @unreachable
end

function Base.:∘(f::CellFieldFromOperation{typeof(get_outward)},ϕ::CellMap)
  @unreachable
end

function Base.:∘(f::CellFieldFromOperation{typeof(get_inward)},ϕ::SkeletonFaceMap)
  left, = merge_cell_fields_at_skeleton(f.args[1]∘ϕ.left,f.args[1]∘ϕ.right)
  left
end

function Base.:∘(f::CellFieldFromOperation{typeof(get_outward)},ϕ::SkeletonFaceMap)
  _, right = merge_cell_fields_at_skeleton(f.args[1]∘ϕ.left,f.args[1]∘ϕ.right)
  right
end

"""
    jump(f)
"""
function jump(sf)
  sf.⁺ - sf.⁻
end

"""
    mean(f)
"""
function mean(sf)
  operate(_mean,sf.⁺,sf.⁻)
end

_mean(x,y) = 0.5*x + 0.5*y

function merge_cell_dofs_at_skeleton(idsL,idsR,axesL,axesR)
  blocks = (idsL,idsR)
  blockids = [(1,),(2,)]
  axs = create_array_of_blocked_axes(axesL,axesR)
  VectorOfBlockArrayCoo(blocks,blockids,axs)
end

function merge_cell_fields_at_skeleton(cfL,cfR)
  @assert is_basis(cfL) == is_basis(cfR)
  if !is_basis(cfL)
    cfL, cfR
  else
    ax1 = get_cell_axes(cfL)
    ax2 = get_cell_axes(cfR)
    arrL = insert_array_of_bases_in_block(1,get_array(cfL),ax1,ax2)
    cfSL = GenericCellField(arrL,arrL.axes,MetaSizeStyle(cfL))
    arrR = insert_array_of_bases_in_block(2,get_array(cfR),ax1,ax2)
    cfSR = GenericCellField(arrR,arrR.axes,MetaSizeStyle(cfR))
    cfSL, cfSR
  end
end


#"""
#    struct CellFieldFromOperation{S} <: CellField
#      op
#      args::Tuple
#    end
#
#Struct representing a `CellField` defined as an operation between other `CellField` objects.
#This struct allows to compute, e.g., `(a*b)∘ϕ` as `(a∘ϕ)*(b∘ϕ)`, which is the trick needed to 
#implement operations involving `CellField` objects defined on different geometrical entities.
#"""
#struct CellFieldFromOperation{F,S} <: CellField
#  op
#  args::Tuple
#  cell_axes::AbstractArray
#  metasize::Val{S}
#  memo::Dict
#  function CellFieldFromOperation(
#    op,::F,args::Tuple,cell_axes::AbstractArray,metasize::Val{S}) where {F,S}
#    memo = Dict()
#    new{F,S}(op,args,cell_axes,metasize,memo)
#  end
#end
#
#function CellFieldFromOperation(op,args::Tuple,cell_axes::AbstractArray,metasize::Val)
#  CellFieldFromOperation(op,op,args,cell_axes,metasize)
#end
#
#function get_array(a::CellFieldFromOperation)
#  key = :get_array
#  memo = get_memo(a)
#  if !haskey(memo,key)
#    memo[key] = get_array(a.op(a.args...))
#  end
#  memo[key]
#end
#
#get_cell_axes(a::CellFieldFromOperation) = a.cell_axes
#get_memo(a::CellFieldFromOperation) = a.memo
#MetaSizeStyle(::Type{CellFieldFromOperation{F,S}}) where {F,S} = Val{S}()
#
#for op in (:+ ,:-)
#  @eval begin
#    function compute_gradient(a::CellFieldFromOperation{typeof($op)})
#      $op(map(gradient,a.args)...)
#    end
#  end
#end
#
#function compute_gradient(a::CellFieldFromOperation{typeof(*)})
#  _product_rule(a.args[1],a.args[2])
#end
#
## TODO this is very hacky. The reason is that we cannot call `get_array(::CellFieldComposedWithInverseMap)`
## in general.  We need to fix this. The other option is to compute the right gradient for numbers,
## but inefficient
#function _product_rule(a,b)
#    a*gradient(b) + b*gradient(a)
#end
#
#function _product_rule(a::GenericCellField,b)
#  _product_rule_a(a,get_array(a),b)
#end
#
#function _product_rule(a,b::GenericCellField)
#  _product_rule_b(a,b,get_array(b))
#end
#
#function _product_rule(a::GenericCellField,b::GenericCellField)
#  _product_rule(a,get_array(a),b,get_array(b))
#end
#
#function _product_rule_a(a,_a,b)
#    a*gradient(b) + b*gradient(a)
#end
#
#function _product_rule_a(a,::AbstractArray{<:Number},b)
#    a*gradient(b)
#end
#
#function _product_rule_b(a,b,_b)
#    a*gradient(b) + b*gradient(a)
#end
#
#function _product_rule_b(a,b,::AbstractArray{<:Number})
#    b*gradient(a)
#end
#
#function _product_rule(a,_a,b,_b)
#  a*gradient(b) + b*gradient(a)
#end
#
#function _product_rule(a,::AbstractArray{<:Number},b,_b)
#    a*gradient(b)
#end
#
#function _product_rule(a,_a,b,::AbstractArray{<:Number})
#    b*gradient(a)
#end
#
#function _product_rule(a,::AbstractArray{<:Number},b,::AbstractArray{<:Number})
#  @notimplemented
#end
#
#function Base.:∘(f::CellFieldFromOperation,ϕ::CellMap)
#  _compose_cell_field_from_op(f,ϕ)
#end
#
#function _compose_cell_field_from_op(f,ϕ)
#  f.op(map(i->i∘ϕ,f.args)...)
#end
#
#function operate(op,a::CellField)
#  key = objectid(op)
#  memo = get_memo(a)
#  if !haskey(memo,key)
#    memo[key] = _operate_wrapper(op,a)
#  end
#  memo[key]
#end
#
#function operate(op,a::CellField...)
#  _operate_wrapper(op,a...)
#end
#
#function operate(op,cf1::CellField,object)
#  cf2 = convert_to_cell_field(object,length(cf1))
#  operate(op,cf1,cf2)
#end
#
#function operate(op,object,cf2::CellField)
#  cf1 = convert_to_cell_field(object,length(cf2))
#  operate(op,cf1,cf2)
#end
#
#function operate(op,a1::CellField,args...)
#  cfs = map(obj->convert_to_cell_field(obj,length(a1)),args)
#  operate(op,a1,cfs...)
#end
#
#function _operate_wrapper(op,args...)
# #if all( length(first(args)) .== map(length,args) )
# if all( length(get_cell_axes(first(args))) .== map(i->length(get_cell_axes(i)),args) )
#   cell_axes, metasize = _compute_metadata_from_op(op,args...)
# else # TODO. This is hacky.
#    # This can be fixed if f∘inverse_map(ϕ) returns an ExtendedVector
#    # for the case that triggers this branch.
#    cell_axes = get_cell_axes(first(args))
#    metasize = MetaSizeStyle(first(args))
#  end
#  f(args...) = _operate(op,args...)
#  CellFieldFromOperation(f,op,args,cell_axes,metasize)
#end
#
#function _operate(op,args::CellField...)
#  a1 = first(args)
#  @assert all( map( a->length(a) == length(a1), args) )
#  arrs = map(get_array,args)
#  m = operate_arrays_of_fields(Fields._UnimplementedField,op,arrs...)
#  cell_axes, metasize = _compute_metadata_from_op(op,args...)
#  GenericCellField(m,cell_axes,metasize)
#end
#
#function _compute_metadata_from_op(op,args::CellField...)
#  axs = apply(field_operation_axes,map(get_cell_axes,args)...)
#  S = field_operation_metasize(map(get_metasize,args)...)
#  axs, Val(S)
#end
#
#function _compute_metadata_from_op(op,arg::CellField)
#  get_cell_axes(arg), MetaSizeStyle(arg)
#end
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
## Bases-related
#
#function trialize_cell_basis(test::CellField)
#  key = :trialize_cell_basis
#  memo = get_memo(test)
#  if !haskey(memo,key)
#    memo[key] = _trialize_cell_basis(test,MetaSizeStyle(test))
#  end
#  memo[key]
#end
#
#function _trialize_cell_basis(test,metasize)
#  @unreachable
#end
#
#function _trialize_cell_basis(trial,metasize::Val{(1,:)})
#  trial
#end
#
#function _trialize_cell_basis(test::CellField,metasize::Val{(:,)})
#  axs = apply(Fields._add_singleton_block,get_cell_axes(test))
#  metasize = Val((1,:))
#  CellFieldFromOperation(trialize_cell_basis,(test,),axs,metasize) do test
#    axs = apply(Fields._add_singleton_block,get_cell_axes(test))
#    metasize = Val((1,:))
#    array = trialize_array_of_bases(get_array(test))
#    GenericCellField(array,axs,metasize)
#  end
#end
#
#function compute_gradient(cf::CellFieldFromOperation{typeof(trialize_cell_basis)})
#  trialize_cell_basis(gradient(first(cf.args)))
#end
#
#function Fields.lincomb(a::CellField,b::AbstractArray)
#  cell_axes = Fill((),length(a))
#  metasize = Val(())
#  CellFieldFromOperation(lincomb,(a,),cell_axes,metasize) do a
#    array = lincomb(get_array(a),b)
#    GenericCellField(array)
#  end
#end
#
#function Fields.lincomb(cf::CellFieldFromOperation{typeof(trialize_cell_basis)},b::AbstractArray)
#  lincomb(first(cf.args),b)
#end
#
#function compute_gradient(cf::CellFieldFromOperation{typeof(lincomb)})
#  a = gradient(cf.args[1])
#  CellFieldFromOperation(cf.op,lincomb,(a,),get_cell_axes(cf),MetaSizeStyle(cf))
#end
#
#function _compose_cell_field_from_op(f::CellFieldFromOperation{typeof(lincomb)},ϕ)
#  a = f.op(f.args...)
#  a∘ϕ
#end
#
#
##function _to_ref_face_space(f,ϕ)
##  cell_to_f = f∘ϕ.cell_map
##  face_to_f = reindex(cell_to_f,ϕ.face_to_cell)
##  face_to_f∘ϕ.refface_to_refcell_map
##end
#
##function Base.:∘(f::CellFieldComposedWithInverseMap,ϕ::ReindexedCellMap)
##  if f.map.direct_cell_map == ϕ
##    return f.f
##  elseif f.map.direct_cell_map == ϕ.cell_map
##    return reindex(f.f,ϕ.ids)
##  else
##    return (f.f∘get_inverse_cell_map(f.map))∘ϕ
##  end
##end
##
##
##function Base.:∘(f::CellFieldComposedWithInverseMap,ϕ::AppendedCellField)
##  if f.map.direct_cell_map == ϕ
##    return f.f
##  else
##    return lazy_append(f∘ϕ.a,f∘ϕ.b)
##  end
##end
##
#
#
#
#
### Lazy append
##
##function Arrays.lazy_append(a::CellField,b::CellField)
##  AppendedCellField(a,b)
##end
##
##struct AppendedCellField{F<:CellField} <: CellField
##  a::F
##  b::CellField
##  memo::Dict
##  function AppendedCellField(a::CellField,b::CellField)
##    F = typeof(a)
##    new{F}(a,b,Dict())
##  end
##end
##
##Arrays.get_array(f::AppendedCellField) = lazy_append(get_array(f.a),get_array(f.b))
##get_cell_axes(f::AppendedCellField) = lazy_append(get_cell_axes(f.a),get_cell_axes(f.b))
##get_memo(f::AppendedCellField) = f.memo
##MetaSizeStyle(::Type{AppendedCellField{T}}) where T = MetaSizeStyle(T)
##Base.:∘(f::CellFieldFromOperation,ϕ::AppendedCellField) = _compose_cell_field_from_op(f,ϕ)
##Base.:∘(f::CellField,ϕ::AppendedCellField) = lazy_append(f∘ϕ.a,f∘ϕ.b)
##Base.:∘(f,ϕ::AppendedCellField) = lazy_append(f∘ϕ.a,f∘ϕ.b)
##
### Other relevant types of maps
##
##struct ReindexedCellMap{F<:CellField} <: CellField
##  map::F
##  cell_map::CellField
##  ids::AbstractArray
##end
##
##function ReindexedCellMap(cell_map::CellField,ids::AbstractArray)
##  map = reindex(cell_map,ids)
##  ReindexedCellMap(map,cell_map,ids)
##end
##
##Arrays.get_array(a::ReindexedCellMap) = get_array(a.map)
##get_cell_axes(a::ReindexedCellMap) = get_cell_axes(a.map)
##get_memo(a::ReindexedCellMap) = get_memo(a.map)
##MetaSizeStyle(::Type{ReindexedCellMap{T}}) where T = MetaSizeStyle(T)
##Base.:∘(f::CellFieldFromOperation,ϕ::ReindexedCellMap) = _compose_cell_field_from_op(f,ϕ)
##
##function Base.:∘(f::CellField,ϕ::ReindexedCellMap)
##  reindex(f∘ϕ.cell_map,ϕ.ids)
##end
##
##struct FaceMap{T<:CellField} <: CellField
##  face_map::T
##  cell_map::CellField
##  face_to_cell::AbstractArray
##  refface_to_refcell_map::CellField
##end
##
##Arrays.get_array(a::FaceMap) = get_array(a.face_map)
##get_cell_axes(a::FaceMap) = get_cell_axes(a.face_map)
##get_memo(a::FaceMap) = get_memo(a.face_map)
##MetaSizeStyle(::Type{FaceMap{T}}) where T = MetaSizeStyle(T)
##Base.:∘(f::CellFieldFromOperation,ϕ::FaceMap) = _compose_cell_field_from_op(f,ϕ)
##function Base.:∘(f::CellField,ϕ::FaceMap)
##  g = reindex(f,ϕ.face_to_cell)
##  g∘ϕ.face_map
##end
##


