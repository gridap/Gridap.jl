
"""
    abstract type GridapType end
"""
abstract type GridapType end

function Base.show(io::IO,object::GridapType)
  print(io,"$(nameof(typeof(object)))()")
end

"""
"""
function unary_operation(op,object::GridapType)
  @unreachable "Unary operation $op is not defined for Gridap objects of type $(typeof(object))"
end

"""
"""
function get_binary_operation(object::GridapType)
  @unreachable "Binary operation rule is not defined for Gridap objects of type $(typeof(object))"
end

function convert_to_operable(a::GridapType,b)
  @unreachable "Not defined how to converto an object of type $(typeof(b)) to another one that is operable with objects of type $(typeof(a))"
end

for op in (:+,:-,:tr, :transpose, :adjoint)
  @eval begin
    function ($op)(a::GridapType)
      unary_operation($op,a)
    end
  end
end

for op in (:+,:-,:*,:cross,:dot)
  @eval begin

    function ($op)(a::GridapType,b::GridapType)
      binop = get_binary_operation(a)
      binop($op,a,b)
    end
    
    function ($op)(a,b::GridapType)
      binop = get_binary_operation(b)
      _a = convert_to_operable(b,a)
      binop($op,_a,b)
    end
    
    function ($op)(a::GridapType,b)
      binop = get_binary_operation(a)
      _b = convert_to_operable(a,b)
      binop($op,a,_b)
    end

  end
end
