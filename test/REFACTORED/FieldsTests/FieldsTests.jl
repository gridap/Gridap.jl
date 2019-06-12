module FieldsTests

using Test
using Gridap

using ..FieldsMocks
using ..FieldsMocks: GradMockField

f = MockField(2,Int)

g = ∇(f)

@test isa(g,GradMockField)

end # module
