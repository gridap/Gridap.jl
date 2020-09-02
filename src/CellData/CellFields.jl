
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

# Tester

"""
    test_cell_field(
      cf::CellField,
      x::AbstractArray,
      b::AbstractArray,
      pred=(==);
      grad=nothing)
"""
function test_cell_field(cf::CellField,x::AbstractArray,b::AbstractArray,pred=(==);grad=nothing)
  a = evaluate(cf,x)
  test_array(a,b,pred)
  if grad != nothing
    g = evaluate(gradient(cf),x)
    test_array(g,grad,pred)
  end
  @test get_metasize(_cf) == get_metasize(cf)
  @test get_cell_axes(_cf) == get_cell_axes(cf)
  @test isa(get_cell_axes(cf),AbstractArray{<:Tuple})
  @test isa(get_metasize(cf),Tuple)
  @test isa(MetaSizeStyle(cf),Val)
  @test isa(get_memo(cf),Dict)
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

#function Base.:∘(a::GenericCellField,ϕ::CellField)
#  _compose_cell_fields(a,get_array(a),ϕ)
#end
#
#function _compose_cell_fields(a::GenericCellField,array,ϕ)
#  _compose_cell_fields(a,ϕ)
#end
#
#function _compose_cell_fields(a::GenericCellField,::AbstractArray{<:Number},ϕ)
#  @assert length(a) == length(ϕ)
#  a
#end
#
#function _compose_cell_fields(a::GenericCellField,f::Fill{<:Fields.FunctionField},ϕ)
#  array = compose(f.value.f,get_array(ϕ))
#  GenericCellField(array)
#end

# The rest of the file is some default API using the CellField interface

"""
    Base.:(==)(a::CellField,b::CellField)

Check if two `CellField` objects are the same. Uses `===` by default. 
"""
function Base.:(==)(a::CellField,b::CellField)
  a === b || get_array(a) === get_array(b)
end

"""
    length(cf::CellField)
"""
function Base.length(cf::CellField)
  length(get_cell_axes(cf))
end

"""
    evaluate(f::CellField,x::AbstractArray)

It is assumed that `x` is in the domain of the array of fields
in `get_array(f)`.
"""
function evaluate(f::CellField,x::AbstractArray)
  key = (:evaluate,objectid(x))
  memo = get_memo(cf)
  if !haskey(memo,key)
    memo[key] = evaluate_array_of_fields(get_array(f),x)
  end
  memo[key]
end

"""
    (f::CellField)(x)

Functor like-evaluation `f(x)` for `CellField` objects. Syntactic sugar for `evaluate(f,x)`.
"""
function (f::CellField)(x)
  evaluate(f,x)
end

"""
    Base.:∘(f::CellField,ϕ::CellField)

Composition of cell fields `(f∘ϕ)(q) == f(ϕ(q)) `
"""
function Base.:∘(f::CellField,ϕ::CellField)
  _compose_cell_fields(f,ϕ)
end

function _compose_cell_fields(f,ϕ)
  @notimplementedif length(f) == length(g)
  array = compose_field_arrays(get_array(f),get_array(ϕ))
  GenericCellField(array,get_cell_axes(f),MetaSizeStyle(f))
end

function Base.:∘(object,ϕ::CellField)
  f = convert_to_cell_field(object,length(ϕ))
  f∘ϕ
end

"""
    reindex(f::CellField,a::AbstractVector)
"""
function Arrays.reindex(f::CellField,a::AbstractVector)
  array = reindex(get_array(f),a)
  cell_axes = reindex(get_cell_axes(f),a)
  GenericCellField(array,cell_axes,MetaSizeStyle(f))
end

function Arrays.reindex(f::object,a::AbstractVector)
  f = convert_to_cell_field(object,length(a))
  reindex(f,a)
end

"""
    gradient(cf::CellField)
"""
function gradient(cf::CellField)
  key = :gradient
  memo = get_memo(cf)
  if !haskey(memo,key)
    a = get_array(cf)
    ag = field_array_gradient(a)
    memo[key] = GenericCellField(ag,get_cell_axes(cf),MetaSizeStyle(cf))
  end
  memo[key]
end

# Operations

"""
    struct CellFieldFromOperation{S} <: CellField
      op
      args::Tuple
    end

Struct representing a `CellField` defined as an operation between other `CellField` objects.
This struct allows to compute, e.g., `(a*b)∘ϕ` as `(a∘ϕ)*(b∘ϕ)`, which is the trick needed to 
implement operations involving `CellField` objects defined on different geometrical objects.
"""
struct CellFieldFromOperation{F,S} <: CellField
  op::F
  args::Tuple
  cell_axes::AbstractArray
  metasize::Val{S}
  memo::Dict
  function CellFieldFromOperation(
    op::F,args::Tuple,cell_axes::AbstractArray,metasize::Val{S}) where {F,S}
    memo = Dict()
    new{F,S}(op,args,cell_axes,metasize,memo)
  end
end

function get_array(a::CellFieldFromOperation)
  key = :get_array
  memo = get_memo(cf)
  if !haskey(memo,key)
    memo[key] = get_array(a.op(a.args...))
  end
  memo[key]
end

get_cell_axes(a::CellFieldFromOperation) = a.cell_axes
get_memo(a::CellFieldFromOperation) = a.memo
MetaSizeStyle(::Type{CellFieldFromOperation{F,S}}) where {F,S} = Val{S}()

function Base.:∘(f::CellFieldFromOperation,ϕ::CellField)
  f.op(map(i->i∘ϕ,f.args)...)
end

function operate(op,a::CellField)
  key = objectid(op)
  memo = get_memo(cf)
  if !haskey(memo,key)
    memo[key] = _operate_wrapper(op,a)
  end
  memo[key]
end

function operate(op,a::CellField...)
  _operate_wrapper(op,a...)
end

function operate(op,cf1::CellField,object)
  cf2 = convert_to_cell_field(object,length(cf1))
  operate(op,cf1,cf2)
end

function operate(op,object,cf1::CellField)
  cf1 = convert_to_cell_field(object,length(cf2))
  operate(op,cf1,cf2)
end

function operate(op,a1::CellField,args...)
  cfs = map(obj->convert_to_cell_field(obj,length(a1)),args)
  operate(op,a1,cfs...)
end

struct Operate{F}
  op::F
end

(a::Operate)(args...) = _operate(a.op,args...)

function _operate_wrapper(op,args...)
  cell_axes, metasize = _compute_metadata_from_op(op,args...)
  CellFieldFromOperation(Operate(op),args,cell_axes,metasize)
end

function _operate(op,args::CellField...)
  a1 = first(args)
  @assert all( map( a->length(a) == length(a1), args) )
  arrs = map(get_array,args)
  m = operate_arrays_of_fields(Fields._UnimplementedField,op,arrs...)
  cell_axes, metasize = _compute_metadata_from_op(op,args...)
  GenericCellField(m,cell_axes,metasize)
end

function _compute_metadata_from_op(op,args::CellField...)
  axs = apply(field_operation_axes,map(get_cell_axes,args)...)
  S = field_operation_metasize(map(get_metasize,args)...)
  axs, Val(S)
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
  axs = apply(Fields._add_singleton_block,get_cell_axes(test))
  metasize = Val((1,:))
  CellFieldFromOperation((test,),axs,metasize) do test
    axs = apply(Fields._add_singleton_block,get_cell_axes(test))
    metasize = Val((1,:))
    array = trialize_array_of_bases(get_array(test))
    GenericCellField(array,axs,metasize)
  end
end

function Fields.lincomb(a::CellField,b::AbstractArray)
  @notimplementedif !is_test(a)
  array = lincomb(get_array(a),b)
  GenericCellField(array)
end

# Skeleton related

"""
    jump(sf)
"""
function jump(sf)
  sf.⁺ - sf.⁻
end

"""
    mean(sf::SkeletonCellField)
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
    cfSL = similar_object(cfL,arrL,arrL.axes,MetaSizeStyle(cfL))
    arrR = insert_array_of_bases_in_block(2,get_array(cfR),ax1,ax2)
    cfSR = similar_object(cfR,arrR,arrR.axes,MetaSizeStyle(cfR))
    cfSL, cfSR
  end
end

# Optimizations associated with inverse maps (i.e., with the reference space)

"""
    inverse_map(ϕ::CellField)

Return another `CellField` representing the inverse map of the given map `ϕ`.
By defalut, it returns an instance of `InverseCellMap`, which keeps the direct cell map `ϕ` as metadata.
"""
inverse_map(ϕ::CellField) = InverseCellMap(ϕ)

inverse_map(ϕ::AbstractArray) = @notimplemented

"""
    struct InverseCellMap <: CellField
      direct_cell_map::CellField
      # + private fields
    end

An inverse cell map that keeps the direct cell map as metadata.
"""
struct InverseCellMap <: CellField
  direct_cell_map::CellField
  memo::Dict
  function InverseCellMap(direct_cell_map::CellField)
    memo = Dict()
    new(direct_cell_map,memo)
  end
end

function get_inverse_cell_map(a::InverseCellMap)
  key = :inverse_cell_map
  memo = get_memo(a)
  if !haskey(memo,key)
    memo[key] = _get_inverse_map(a)
  end
  memo[key]
end

function _get_inverse_map(a)
  array = inverse_map(get_array(a.direct_cell_map))
  @notimplementedif length(array) != length(direct_cell_map)
  GenericCellField(
    array,get_cell_axes(a.direct_cell_map),MetaSizeStyle(a.direct_cell_map))
end

get_array(a::InverseCellMap) = get_array(get_inverse_cell_map(a))
get_cell_axes(a::InverseCellMap) = get_cell_axes(a.direct_cell_map)
get_memo(a::InverseCellMap) = get_memo(a.memo)
MetaSizeStyle(::Type{InverseCellMap}) = Val(())

function Arrays.reindex(f::InverseCellMap,a::AbstractVector)
  InverseCellMap(reindex(f.direct_cell_map,a),reindex(f.inverse_cell_map))
end

"""
    struct CellFieldComposedWithInverseMap{T<:CellField} <: CellField
      f::T
      map::InverseCellMap
    end
struct representing `f∘inverse_map(ϕ)`, i.e. a `CellField` in the reference space
pushed to the physical one.  This struct allows a number of optimizations.
In particular, the computation of `f∘inverse_map(ϕ)∘ϕ`, which is simply computed as`f`.
It also computes `∇(f∘inverse_map(ϕ))` as `(inv(ϕ)*∇(f))∘inverse_map(ϕ)` which is in turn
stored in an instance of `CellFieldComposedWithInverseMap`.
"""
struct CellFieldComposedWithInverseMap{T<:CellField} <: CellField
  f::T
  map::InverseCellMap
end

function get_array(a::CellFieldComposedWithInverseMap)
  get_array(f.f∘f.map.inverse_map)
end

get_cell_axes(a::CellFieldComposedWithInverseMap) = @notimplemented
get_memo(a::CellFieldComposedWithInverseMap) = @notimplemented
MetaSizeStyle(::Type{CellFieldComposedWithInverseMap}) = @notimplemented
evaluate(a::CellFieldComposedWithInverseMap,x::AbstractArray) = @notimplemented

function Base.:∘(f::CellFieldComposedWithInverseMap,ϕ::CellField)
  if f.map.direct_cell_map == ϕ
    return f.f
  else
    @notimplemented
    return (f.f∘f.map.inverse_map)∘ϕ
  end
end

function Arrays.reindex(f::CellFieldComposedWithInverseMap,a::AbstractVector)
  CellFieldComposedWithInverseMap(reindex(f.f,a),reindex(f.map,a))
end

# TODO
function gradient(cf::CellFieldComposedWithInverseMap)
  @notimplemented
end

function Fields.lincomb(a::CellFieldComposedWithInverseMap,b::AbstractArray)
  u = lincomb(a.f,b)
  CellFieldComposedWithInverseMap(u,a.map)
end
