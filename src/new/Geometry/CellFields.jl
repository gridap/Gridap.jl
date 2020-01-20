
"""
"""
abstract type CellFieldLike <: GridapType end

"""
"""
function get_array(cf::CellFieldLike)
  @abstractmethod
end

"""
"""
function get_cell_map(cf::CellFieldLike)
  @abstractmethod
end

"""
"""
function gradient(cf::CellFieldLike)
  @abstractmethod
end

"""
"""
function grad2curl(cf::CellFieldLike)
  @abstractmethod
end

"""
"""
function test_cell_field_like(cf::CellFieldLike,x::AbstractArray,b::AbstractArray,pred=(==);grad=nothing)
  cell_map = get_cell_map(cf)
  @test isa(cell_map,AbstractArray)
  a = evaluate(cf,x)
  test_array(a,b,pred)
  if grad != nothing
    g = evaluate(gradient(cf),x)
    test_array(g,grad,pred)
  end
end

"""
"""
function evaluate(cf::CellFieldLike,x)
  a = get_array(cf)
  evaluate_field_array(a,x)
end

"""
"""
function Base.length(cf::CellFieldLike)
  a = get_array(cf)
  length(a)
end

"""
"""
abstract type CellField <: CellFieldLike end

"""
"""
function test_cell_field(cf::CellField,args...;kwargs...)
  test_cell_field_like(cf,args...;kwargs...)
end

"""
"""
function similar_object(cf::CellField,array::AbstractArray)
  cm = get_cell_map(cf)
  GenericCellField(array,cm)
end

function similar_object(cf1::CellField,cf2::CellField,array::AbstractArray)
  cm = get_cell_map(cf1)
  GenericCellField(array,cm)
end

# Diff operations

struct UnimplementedField <: Field end

function gradient(cf::CellField)
  a = get_array(cf)
  g = field_array_gradient(a)
  similar_object(cf,g)
end

function grad2curl(cf::CellField)
  a = get_array(cf)
  g = grad2curl(UnimplementedField,a)
  similar_object(cf,g)
end

# Operations

function operate(op,cf::CellField)
  a = get_array(cf)
  b = field_array_operation(UnimplementedField,op,a)
  similar_object(cf,b)
end

function operate(op,cf1::CellField,cf2::CellField)
  @assert length(cf1) == length(cf2)
  a1 = get_array(cf1)
  a2 = get_array(cf2)
  b = field_array_operation(UnimplementedField,op,a1,a2)
  similar_object(cf1,cf2,b)
end

function operate(op,cf1::CellField,object::Union{Function,Number})
  cm = get_cell_map(cf1)
  cf2 = convert_to_cell_field(object,cm)
  operate(op,cf1,cf2)
end

function operate(op,object::Union{Function,Number},cf2::CellField)
  cm = get_cell_map(cf2)
  cf1 = convert_to_cell_field(object,cm)
  operate(op,cf1,cf2)
end

# Conversions

"""
"""
function convert_to_cell_field(object::CellField,cell_map)
  object
end

function convert_to_cell_field(object::AbstractArray,cell_map)
  @assert length(object) == length(cell_map)
  GenericCellField(object,cell_map)
end

function convert_to_cell_field(object::Function,cell_map)
  b = compose(object,cell_map)
  GenericCellField(b,cell_map)
end

function convert_to_cell_field(object::Number,cell_map)
  a = Fill(object,length(cell_map))
  GenericCellField(a,cell_map)
end

# Concrete implementation

"""
"""
struct GenericCellField{A<:AbstractArray,B<:AbstractArray} <: CellField
  array::A
  cell_map::B
end

function get_array(cf::GenericCellField)
  cf.array
end

function get_cell_map(cf::GenericCellField)
  cf.cell_map
end


# Skeleton related

"""
"""
struct SkeletonCellField <: GridapType
  left::CellField
  right::CellField
end

"""
"""
function get_cell_map(a::SkeletonCellField)
  get_cell_map(a.left)
end

"""
"""
function jump(sf::SkeletonCellField)
  sf.left - sf.right
end

"""
"""
function mean(sf::SkeletonCellField)
  operate(_mean,sf.left,sf.right)
end

_mean(x,y) = 0.5*x + 0.5*y

function gradient(cf::SkeletonCellField)
  left = gradient(cf.left)
  right = gradient(cf.right)
  SkeletonCellField(left,right)
end

function grad2curl(cf::SkeletonCellField)
  left = grad2curl(cf.left)
  right = grad2curl(cf.right)
  SkeletonCellField(left,right)
end

"""
"""
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

function operate(op,cf1::SkeletonCellField,object::Union{Function,Number})
  cm = get_cell_map(cf1)
  cf2 = convert_to_cell_field(object,cm)
  operate(op,cf1,cf2)
end

function operate(op,object::Union{Function,Number},cf2::SkeletonCellField)
  cm = get_cell_map(cf2)
  cf1 = convert_to_cell_field(object,cm)
  operate(op,cf1,cf2)
end

