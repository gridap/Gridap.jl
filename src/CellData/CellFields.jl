
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
    get_cell_map(cf::CellField)
"""
function get_cell_map(cf::CellField)
  @abstractmethod
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
This trait returns `Val{true}()` when the `CellField` is defined in a
reference finite element space, and `Val{false}()` when it is defined in the
physical space
"""
RefStyle(::Type{<:CellField}) = @notimplemented

# We use duck typing here for all types marked with the RefStyle
RefStyle(::Type) = @notimplemented
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
  cell_map = get_cell_map(cf)
  @test isa(cell_map,AbstractArray)
  a = evaluate(cf,x)
  test_array(a,b,pred)
  if grad != nothing
    g = evaluate(gradient(cf),x)
    test_array(g,grad,pred)
  end
  rs = RefStyle(cf)
  @test isa(get_val_parameter(rs),Bool)
  _cf = change_ref_style(cf)
  @test get_array(_cf) === get_array(cf)
  @test is_in_ref_space(cf) == !is_in_ref_space(_cf)
  @test is_in_physical_space(cf) == !is_in_physical_space(_cf)
  @test get_metasize(_cf) == get_metasize(cf)
  @test get_cell_axes(_cf) == get_cell_axes(cf)
  @test isa(get_cell_axes(cf),AbstractArray{<:Tuple})
  @test isa(get_metasize(cf),Tuple)
  @test isa(MetaSizeStyle(cf),Val)
  @test isa(get_memo(cf),Dict)
end

# Concrete implementation

"""
struct GenericCellField{R,S} <: CellField
  # Private fields
end

"""
struct GenericCellField{R,S} <: CellField
  array::AbstractArray
  cell_map::AbstractArray
  ref_trait::Val{R}
  cell_axes::AbstractArray
  metasize::Val{S}
  memo::Dict
  function GenericCellField(
    array::AbstractArray,
    cell_map::AbstractArray,
    ref_trait::Val{R},
    cell_axes::AbstractArray,
    metasize::Val{S}) where {R,S}
    memo = Dict()
    new{R,S}(array,cell_map,ref_trait,cell_axes,metasize,memo)
  end
end

function GenericCellField(array::AbstractArray,cell_map::AbstractArray)
  GenericCellField(array,cell_map,Val{true}(),Fill((),length(cell_map)))
end

function GenericCellField(array::AbstractArray,cell_map::AbstractArray,ref_trait::Val)
  GenericCellField(array,cell_map,ref_trait,Fill((),length(cell_map)))
end

function GenericCellField(
  array::AbstractArray,cell_map::AbstractArray,ref_trait::Val,cell_axes::AbstractArray)
  GenericCellField(array,cell_map,ref_trait,cell_axes,Val(()))
end

function get_array(cf::GenericCellField)
  cf.array
end

function get_cell_map(cf::GenericCellField)
  cf.cell_map
end

function get_memo(cf::GenericCellField)
  cf.memo
end

function get_cell_axes(cf::GenericCellField)
  cf.cell_axes
end

function RefStyle(::Type{<:GenericCellField{R}}) where {R}
  Val{R}()
end

function MetaSizeStyle(::Type{<:GenericCellField{R,S}}) where {R,S}
  Val{S}()
end

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
  @assert is_test(test)
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

