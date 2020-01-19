module GridapTypesTests

using Test
using Gridap.Helpers

struct Foo{A} <: GridapType
  bar::A
end

object = Foo(nothing)

@test sprint(show,"text/plain",object) == "$(nameof(typeof(object)))()"

end # module
