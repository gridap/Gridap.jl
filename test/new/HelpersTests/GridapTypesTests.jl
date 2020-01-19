module GridapTypesTests

using Test
using Gridap.Helpers

import Gridap.Helpers: unary_operation
import Gridap.Helpers: get_binary_operation
import Gridap.Helpers: convert_to_operable

struct Foo{A} <: GridapType
  bar::A
end

object = Foo(nothing)

@test sprint(show,"text/plain",object) == "$(nameof(typeof(object)))()"

function unary_operation(op,a::Foo)
  bar = op(a.bar)
  Foo(bar)
end

function get_binary_operation(a::Foo)
  (op,a,b) -> Foo(op(a.bar,b.bar))
end

function convert_to_operable(::Foo,b::Number)
  Foo(b)
end

a = Foo(1)

b = -a
@test b.bar == -1

c = a + b

@test c.bar == 0

c = a - 2
@test c.bar == -1

c = 2 - a
@test c.bar == 1

using Gridap.TensorValues

c = inner(2,a)
@test c.bar == 2

end # module
