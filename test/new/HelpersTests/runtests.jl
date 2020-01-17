module HelpersTests

using Test
using Gridap.Helpers

for D = 0:5
  @test tfill(2,Val(D)) == tuple(fill(2,D)...)
end

struct Foo{A} <: GridapType
  bar::A
end

object = Foo(nothing)

@test sprint(show,"text/plain",object) == "$(nameof(typeof(object)))()"

end # module

