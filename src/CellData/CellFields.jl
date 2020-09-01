
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
struct GenericCellField{R,S} <: CellField
  array::AbstractArray
  cell_axes::AbstractArray
  metasize::Val{S}
  memo::Dict
  function GenericCellField(
    array::AbstractArray,
    cell_axes::AbstractArray,
    metasize::Val{S}) where {S}
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

"""
    Base.:(==)(a::CellField,b::CellField)

Check if two `CellField` objects are the same. Uses `===` by default. 
"""
function Base.:(==)(a::CellField,b::CellField)
  a === b
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
  @notimplementedif length(f) == length(g)
  array = compose_field_arrays(get_array(f),get_array(ϕ))
  GenericCellField(array,get_cell_axes(f),MetaSizeStyle(f))
end

"""
    reindex(f::CellField,a::AbstractVector)
"""
function Arrays.reindex(f::CellField,a::AbstractVector)
  array = reindex(get_array(f),a)
  cell_axes = reindex(get_cell_axes(f),a)
  GenericCellField(array,cell_axes,MetaSizeStyle(f))
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
  inverse_cell_map::CellField
  function InverseCellMap(direct_cell_map::CellField)
    array = inverse_map(get_array(direct_cell_map))
    @notimplementedif length(array) != length(direct_cell_map)
    inverse_cell_map = GenericCellField(
      array,get_cell_axes(direct_cell_map),MetaSizeStyle(direct_cell_map))
    new(direct_cell_map,inverse_cell_map)
  end
end

get_array(a::InverseCellMap) = get_array(a.inverse_cell_map)
get_cell_axes(a::InverseCellMap) = get_cell_axes(a.inverse_cell_map)
get_memo(a::InverseCellMap) = get_memo(a.inverse_cell_map)
MetaSizeStyle(::Type{InverseCellMap}) = Val(())

function Arrays.reindex(f::InverseCellMap,a::AbstractVector)
  InverseCellMap(reindex(f.direct_cell_map,a),reindex(f.inverse_cell_map))
end

"""
    struct CellFieldComposedWithInverseMap{T<:CellField} <: CellField
      f::T
      map::InverseCellMap
    end
struct representing `f∘inverse_map(ϕ)`.  This struct allows a number of optimizations.
In particular, the computation of `f∘inverse_map(ϕ)∘ϕ`, which is simply computed as`f`.
It also computes `∇(f∘inverse_map(ϕ))` as `(inv(ϕ)*∇(f))∘inverse_map(ϕ)` which is in turn
stored in an instance of `CellFieldComposedWithInverseMap`.
"""
struct CellFieldComposedWithInverseMap{T<:CellField} <: CellField
  f::T
  map::InverseCellMap
end

# TODO
get_array(a::CellFieldComposedWithInverseMap) = @notimplemented
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
  CellFieldComposedWithInverseMap(reindex(f.f,a),reindex(a.map))
end

# TODO
function gradient(cf::CellFieldComposedWithInverseMap)
  @notimplemented
end

"""
    struct CellFieldFromOperation{S} <: CellField
      op
      args::Tuple
    end

Struct representing a `CellField` defined as an operation between other `CellField` objects.
This struct allows to compute, e.g., `(a*b)∘ϕ` as `(a∘ϕ)*(b∘ϕ)`, which is the trick needed to 
implement operations involving `CellField` objects defined on different geometrical objects.
"""
struct CellFieldFromOperation{S} <: CellField
  op
  args::Tuple
end

# They are not needed in practice. So, we don't implement them for the moment.
get_array(a::CellFieldFromOperation) = @notimplemented
get_cell_axes(a::CellFieldFromOperation) = @notimplemented
get_memo(a::CellFieldFromOperation) = @notimplemented
MetaSizeStyle(::Type{CellFieldFromOperation}) = @notimplemented
evaluate(a::CellFieldFromOperation,x::AbstractArray) = @notimplemented











get_array(cf,::Val{false}) = @abstractmethod

function get_array(f,::Val{true})
  ϕ = get_cell_map(f)
  get_array((f∘ϕ)∘inverse_map(ϕ))
end

inverse_map(ϕ::CellField) = inverse_map(get_array(ϕ))
inverse_map(ϕ::AbstractArray) = @notimplemented

function Base.:∘(f::CellField,ϕ::CellField)
  _compose(f,RefStyle(f),ϕ)
end

function _compose(f,::Val{false},ϕ)
  a = compose_field_arrays(get_array(f),get_array(ϕ))
  similar_object(f,a) # TODO
end

function _compose(f,::Val{true},ϕ)
  if get_cell_map(f) == ϕ 
    similar_object(f,get_ref_array(f)) # TODO
  else
    _compose(f,Val(false),ϕ)
  end
end

function Base.:(==)(a::CellField,b::CellField)
  a === b
end

"""
    get_cell_axes(cf::CellField)
"""
function get_cell_axes(cf::CellField)
  @abstractmethod
end

function get_memo(cf::CellField)
  @abstractmethod
end

"""
This trait returns `Val(true)` when the concept of reference space
is defined for the `CellField`, and `Val(false)` otherwise.
All `CellField` types such that this trait is `Val(true)` should define a function
[`get_cell_map(cf::CellField)`](@ref) which returns the cell map defining the reference space.
A cell map is a `CellField` object whose domain is the reference space.
Moreover, types with this trait need to implement
`get_ref_array(cf::CellField)` that returns an array of fields defined in the reference space.
In this case, `get_array(cf::CellField)` returns an array of fields in the physical space.
"""
RefStyle(::Type{<:CellField}) = @notimplemented

"""
    get_cell_map(cf::CellField)
"""
function get_cell_map(cf::CellField)
  @abstractmethod
end

function get_ref_array(cf::CellField)
  @abstractmethod
end

# We use duck typing here for all types marked with the RefStyle
RefStyle(::Type) = @abstractmethod
RefStyle(a) = RefStyle(typeof(a))
is_in_ref_space(::Type{T}) where T = get_val_parameter(RefStyle(T))
is_in_ref_space(::T) where T = is_in_ref_space(T)
is_in_physical_space(::Type{T}) where T = !is_in_ref_space(T)
is_in_physical_space(a::T) where T = !is_in_ref_space(T)

"""
This trait provides information about the field size at a single
evaluation point.

For physical fields `MetaSizeStyle(T) == Val( () )`
For test bases `MetaSizeStyle(T) == Val( (:,) )`
For trial bases `MetaSizeStyle(T) = Val( (1,:) )`
For fields representing elemental matrices `MetaSizeStyle(T) = Val( (:,:) )`
"""
MetaSizeStyle(::Type{<:CellField}) = @notimplemented

# We use duck typing here for all types marked with the MetaSizeStyle trait
MetaSizeStyle(::Type) = @notimplemented
MetaSizeStyle(::T) where T = MetaSizeStyle(T)
get_metasize(::Type{T}) where T = get_val_parameter(MetaSizeStyle(T))
get_metasize(::T) where T = get_val_parameter(MetaSizeStyle(T))
is_test(::Type{T}) where T = get_metasize(T) == (:,)
is_test(::T) where T = is_test(T)
is_trial(::Type{T}) where T = get_metasize(T) == (1,:)
is_trial(::T) where T = is_trial(T)
is_basis(::Type{T}) where T =  is_test(T) || is_trial(T)
is_basis(::T) where T = is_basis(T)



# The rest of the file is some default API using the CellField interface

# Preserving and manipulating the metadata

"""
    similar_object(cf::CellField,array::AbstractArray,cell_axes::AbstractArray,msize_style::Val)
"""
function similar_object(cf::CellField,array::AbstractArray,cell_axes::AbstractArray,msize_style::Val)
  cm = get_cell_map(cf)
  similar_object(cf,array,cm,cell_axes,msize_style)
end

function similar_object(
  cf::CellField,array::AbstractArray,cell_map::AbstractArray,cell_axes::AbstractArray,msize_style::Val)
  GenericCellField(array,cell_map,RefStyle(cf),cell_axes,msize_style)
end

"""
   change_ref_style(cf::CellField)

Return an object with the same array and metadata as in `cf`, except for `RefStyle` which is changed.
"""
function change_ref_style(cf::CellField)
  ref_sty = RefStyle(cf)
  bool = !get_val_parameter(ref_sty)
  new_sty = Val{bool}()
  ar = get_array(cf)
  cm = get_cell_map(cf)
  GenericCellField(ar,cm,new_sty,get_cell_axes(cf),MetaSizeStyle(cf))
end

# Move between ref and phys spaces by using the underling cell_map

to_ref_space(a::CellField) = _to_ref_space(a,RefStyle(a))
_to_ref_space(a,::Val{true}) = a
function _to_ref_space(a,::Val{false})
  cell_map = get_cell_map(a)
  array = compose(  get_array(a), cell_map  )
  no = similar_object(a,array,get_cell_axes(a),MetaSizeStyle(a))
  change_ref_style(no)
end

to_physical_space(a::CellField) = _to_physical_space(a,RefStyle(a))
_to_physical_space(a,::Val{true}) = @notimplemented # and probably not doable in some cases
_to_physical_space(a,::Val{false}) = a

# Assumption : x ALWAIS defined in the reference space
# In the future we can also add the RefStyle to x by defining CellPoints
"""
    evaluate(cf::CellField,x)
"""
function evaluate(cf::CellField,x::AbstractArray)
  key = (:evaluate,objectid(x))
  memo = get_memo(cf)
  if !haskey(memo,key)
    memo[key] = _evaluate(cf,x,RefStyle(cf))
  end
  memo[key]
end

function _evaluate(cf::CellField,x::AbstractArray,::Val{true})
  evaluate_field_array(get_array(cf),x)
end

function _evaluate(cf::CellField,x::AbstractArray,::Val{false})
  cm = get_cell_map(cf)
  _x = evaluate(cm,x)
  evaluate_field_array(get_array(cf),_x)
end

"""
    length(cf::CellField)
"""
function Base.length(cf::CellField)
  a = get_array(cf)
  length(a)
end

function Arrays.reindex(cf::CellField,a::AbstractVector)
  array = reindex(get_array(cf),a)
  cell_axes = reindex(get_cell_axes(cf),a)
  cell_map = reindex(get_cell_map(cf),a)
  similar_object(cf,array,cell_map,cell_axes,MetaSizeStyle(cf))
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
  array = trialize_array_of_bases(get_array(test))
  axs = apply(Fields._add_singleton_block,get_cell_axes(test))
  similar_object(test,array,axs,Val((1,:)))
end

function Fields.lincomb(a::CellField,b::AbstractArray)
  @notimplementedif !is_test(a)
  array = lincomb(get_array(a),b)
  axs = Fill((),length(b))
  similar_object(a,array,axs,Val(()))
end

# Diff ops

"""
    gradient(cf::CellField)
"""
function gradient(cf::CellField)
  key = :gradient
  memo = get_memo(cf)
  if !haskey(memo,key)
    a = get_array(cf)
    ag = field_array_gradient(a)
    memo[key] = similar_object(cf,ag,get_cell_axes(cf),MetaSizeStyle(cf))
  end
  memo[key]
end

# Operations

function operate(op,cf::CellField)
  key = objectid(op)
  memo = get_memo(cf)
  if !haskey(memo,key)
    a = get_array(cf)
    b = operate_arrays_of_fields(Fields._UnimplementedField,op,a)
    memo[key] = similar_object(cf,b,get_cell_axes(cf),MetaSizeStyle(cf))
  end
  memo[key]
end

function operate(op,cf1::CellField,cf2::CellField)
  @assert length(cf1) == length(cf2)
  @assert RefStyle(cf1) == RefStyle(cf2)
  a1 = get_array(cf1)
  a2 = get_array(cf2)
  b = operate_arrays_of_fields(Fields._UnimplementedField,op,a1,a2)
  axs = apply(field_operation_axes,get_cell_axes(cf1),get_cell_axes(cf2))
  metasize = field_operation_metasize(get_metasize(cf1),get_metasize(cf2))
  similar_object(cf1,b,axs,Val(metasize))
end

function operate(op,cf1::CellField,object)
  cm = get_cell_map(cf1)
  cf2 = convert_to_cell_field(object,cm,RefStyle(cf1))
  operate(op,cf1,cf2)
end

function operate(op,object,cf2::CellField)
  cm = get_cell_map(cf2)
  cf1 = convert_to_cell_field(object,cm,RefStyle(cf2))
  operate(op,cf1,cf2)
end

function operate(op,args::CellField...)
  a1 = first(args)
  @assert all( map( a->length(a) == length(a1), args) )
  @assert all( map( a->RefStyle(a)==RefStyle(a1), args) )
  arrs = map(get_array,args)
  m = operate_arrays_of_fields(Fields._UnimplementedField,op,arrs...)
  axs = apply(field_operation_axes,map(get_cell_axes,args)...)
  metasize = field_operation_metasize(map(get_metasize,args)...)
  similar_object(a1,m,axs,Val(metasize))
end

function operate(op,a1::CellField,args...)
  cm = get_cell_map(a1)
  rs = RefStyle(cm)
  cfs = map(obj->convert_to_cell_field(obj,cm,rs),args)
  operate(op,a1,cfs...)
end

# Conversions

"""
    convert_to_cell_field(object,cell_map,ref_style)
"""
function convert_to_cell_field(object,cell_map,ref_style::Val)
  @abstractmethod
end

function convert_to_cell_field(object,cell_map)
  ref_style = Val{true}()
  convert_to_cell_field(object,cell_map,ref_style)
end

function convert_to_cell_field(object::CellField,cell_map)
  @assert length(object) == length(cell_map)
  object
end

# pre-defined conversions

function convert_to_cell_field(object::CellField,cell_map,ref_style::Val)
  @assert RefStyle(object) == ref_style
  object
end

function convert_to_cell_field(object::AbstractArray,cell_map,ref_style::Val)
  @assert length(object) == length(cell_map)
  GenericCellField(object,cell_map,ref_style)
end

function convert_to_cell_field(object::Function,cell_map,ref_style::Val{true})
  b = compose(object,cell_map)
  GenericCellField(b,cell_map,Val{true}())
end

function convert_to_cell_field(fun::Function,cell_map,ref_style::Val{false})
  field = function_field(fun)
  cell_field = Fill(field,length(cell_map))
  GenericCellField(cell_field,cell_map,Val{false}())
end

function convert_to_cell_field(object::Number,cell_map,ref_style::Val)
  a = Fill(object,length(cell_map))
  GenericCellField(a,cell_map,ref_style)
end

# Skeleton related

"""
    struct SkeletonCellField <: GridapType
      left::CellField
      right::CellField
    end

Supports the same differential and algebraic operations than [`CellField`](@ref)
"""
struct SkeletonCellField <: GridapType
  left::CellField
  right::CellField
end

get_inward(a::SkeletonCellField) = a.left

get_outward(a::SkeletonCellField) = a.right

function Base.getproperty(x::SkeletonCellField, sym::Symbol)
  if sym in (:inward,:⁺)
    get_inward(x)
  elseif sym in (:outward,:⁻)
    get_outward(x)
  else
    getfield(x, sym)
  end
end

"""
    get_cell_map(a::SkeletonCellField)
"""
function get_cell_map(a::SkeletonCellField)
  get_cell_map(a.left)
end

"""
    jump(sf::SkeletonCellField)
"""
function jump(sf::SkeletonCellField)
  sf.⁺ - sf.⁻
end

"""
    mean(sf::SkeletonCellField)
"""
function mean(sf::SkeletonCellField)
  operate(_mean,sf.⁺,sf.⁻)
end

_mean(x,y) = 0.5*x + 0.5*y

function gradient(cf::SkeletonCellField)
  left = gradient(cf.left)
  right = gradient(cf.right)
  SkeletonCellField(left,right)
end

function operate(op,cf::SkeletonCellField)
  left = operate(op,cf.left)
  right = operate(op,cf.right)
  SkeletonCellField(left,right)
end

function operate(op,cf1::SkeletonCellField,cf2::SkeletonCellField)
  left = operate(op,cf1.left,cf2.left)
  right = operate(op,cf1.right,cf2.right)
  SkeletonCellField(left,right)
end

function operate(op,cf1::SkeletonCellField,cf2::CellField)
  left = operate(op,cf1.left,cf2)
  right = operate(op,cf1.right,cf2)
  SkeletonCellField(left,right)
end

function operate(op,cf1::CellField,cf2::SkeletonCellField)
  left = operate(op,cf1,cf2.left)
  right = operate(op,cf1,cf2.right)
  SkeletonCellField(left,right)
end

function operate(op,cf1::SkeletonCellField,object)
  cm = get_cell_map(cf1)
  cf2 = convert_to_cell_field(object,cm,RefStyle(cf1.left))
  operate(op,cf1,cf2)
end

function operate(op,object,cf2::SkeletonCellField)
  cm = get_cell_map(cf2)
  cf1 = convert_to_cell_field(object,cm,RefStyle(cf2.left))
  operate(op,cf1,cf2)
end

function merge_cell_dofs_at_skeleton(idsL,idsR,axesL,axesR)
  blocks = (idsL,idsR)
  blockids = [(1,),(2,)]
  axs = create_array_of_blocked_axes(axesL,axesR)
  VectorOfBlockArrayCoo(blocks,blockids,axs)
end

function merge_cell_fields_at_skeleton(cfL,cfR)
  @assert is_basis(cfL) == is_basis(cfR)
  if !is_basis(cfL)
    SkeletonCellField(cfL,cfR)
  else
    ax1 = get_cell_axes(cfL)
    ax2 = get_cell_axes(cfR)
    arrL = insert_array_of_bases_in_block(1,get_array(cfL),ax1,ax2)
    cfSL = similar_object(cfL,arrL,arrL.axes,MetaSizeStyle(cfL))
    arrR = insert_array_of_bases_in_block(2,get_array(cfR),ax1,ax2)
    cfSR = similar_object(cfR,arrR,arrR.axes,MetaSizeStyle(cfR))
    SkeletonCellField(cfSL,cfSR)
  end
end
