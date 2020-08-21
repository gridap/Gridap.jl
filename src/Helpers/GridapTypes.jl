
"""
    abstract type GridapType end
"""
abstract type GridapType end

function Base.show(io::IO,object::GridapType)
  print(io,"$(nameof(typeof(object)))()")
end

"""
"""
function operate(op,a)
  s = "Unary operation $op is not defined for objects of type $(typeof(a))\n"
  s *= "Possible fix, define:\n"
  s *= "Gridap.Helpers.operate(::$(typeof(op)),::$(typeof(a)))"
  @unreachable s
end

"""
"""
function operate(op,a,b)
  s = "Binary operation $op is not defined for objects of type $(typeof(a)) and $(typeof(b))\n"
  s *= "Possible fix, define:\n"
  s *= "Gridap.Helpers.operate(::$(typeof(op)),::$(typeof(a)),::$(typeof(b)))"
  @unreachable s
end

for op in (:+,:-,:tr, :transpose, :adjoint)
  @eval begin
    function ($op)(a::GridapType)
      operate($op,a)
    end
  end
end

for op in (:+,:-,:*,:cross,:dot,:/)
  @eval begin

    function ($op)(a::GridapType,b::GridapType)
      operate($op,a,b)
    end

    function ($op)(a::GridapType,b)
      operate($op,a,b)
    end

    function ($op)(a,b::GridapType)
      operate($op,a,b)
    end

  end
end

