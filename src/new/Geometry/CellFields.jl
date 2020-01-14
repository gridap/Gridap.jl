
"""
"""
abstract type CellField <: GridapType end

"""
"""
function get_array(cf::CellField)
  @abstractmethod
end

"""
"""
function get_cell_map(cf::CellField)
  @abstractmethod
end

"""
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
end

"""
"""
function evaluate(cf::CellField,x)
  a = get_array(cf)
  evaluate_field_array(a,x)
end

"""
"""
function Base.length(cf::CellField)
  a = get_array(cf)
  length(a)
end

"""
"""
function similar_cell_field(cf::CellField,array::AbstractArray)
  cm = get_cell_map(cf)
  GenericCellField(array,cm)
end

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

# Operations

function gradient(cf::CellField)
  a = get_array(cf)
  g = field_array_gradient(a)
  similar_cell_field(cf,g)
end

function grad2curl(cf::CellField)
  a = get_array(cf)
  g = grad2curl(a)
  similar_cell_field(cf,g)
end

(*)(::typeof(∇),f::CellField) = divergence(f)

outer(::typeof(∇),f::CellField) = gradient(f)

outer(f::CellField,::typeof(∇)) = transpose(gradient(f))

cross(::typeof(∇),f::CellField) = curl(f)

struct UnimplementedField <: Field end

"""
"""
function operate_cell_field(op,cf::CellField)
  a = get_array(cf)
  b = field_array_operation(UnimplementedField,op,a)
  similar_cell_field(cf,b)
end

function operate_cell_field(op,cf1::CellField,cf2::CellField)
  @assert length(cf1) == length(cf2)
  a1 = get_array(cf1)
  a2 = get_array(cf2)
  b = field_array_operation(UnimplementedField,op,a1,a2)
  similar_cell_field(cf1,b)
end

for op in (:+,:-,:tr, :transpose, :adjoint, :symmetic_part)
  @eval begin
    function ($op)(cf::CellField)
      operate_cell_field($op,cf)
    end
  end
end

for op in (:+,:-,:*,:inner,:outer)
  @eval begin

    function ($op)(cf1::CellField,cf2::CellField)
      operate_cell_field($op,cf1,cf2)
    end

    function ($op)(cf1::CellField,object)
      cm = get_cell_map(cf1)
      cf2 = convert_to_cell_field(object,cm)
      $op(cf1,cf2)
    end

    function ($op)(object,cf2::CellField)
      cm = get_cell_map(cf2)
      cf1 = convert_to_cell_field(object,cm)
      $op(cf1,cf2)
    end

  end
end

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

