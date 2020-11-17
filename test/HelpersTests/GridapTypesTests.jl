module GridapTypesTests

using Test
using Gridap.Helpers

#import Gridap.Helpers: operate

struct Foo{A} <: GridapType
  bar::A
end

object = Foo(nothing)

@test sprint(show,"text/plain",object) == "$(nameof(typeof(object)))()"

#function operate(op,a::Foo)
#  Foo(op(a.bar))
#end
#
#function operate(op,a::Foo,b::Foo)
#  Foo(op(a.bar,b.bar))
#end
#
#function operate(op,a::Number,b::Foo)
#  operate(op,Foo(a),b)
#end
#
#function operate(op,a::Foo,b::Number)
#  operate(op,a,Foo(b))
#end
#
#a = Foo(1)
#
#b = -a
#@test b.bar == -1
#
#c = a + b
#
#@test c.bar == 0
#
#c = a - 2
#@test c.bar == -1
#
#c = 2 - a
#@test c.bar == 1
#
#using Gridap.TensorValues
#
#c = inner(2,a)
#@test c.bar == 2
#
#c = a / b
#@test c.bar == -1

end # module
