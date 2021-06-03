module InverseFieldsTests

using Gridap.Fields
using Gridap.Arrays
using Gridap.TensorValues

using Test

b0 = Point(0,0)
m0 = TensorValue(1,0,0,1)
id = AffineMap(m0,b0)

b1 = Point(1,1)
m1 = TensorValue(2,0,0,2)
h1 = AffineMap(m1,b1)

b2 = Point(3,3)
m2 = TensorValue(4,0,0,4)
h2 = AffineMap(m2,b2)

h = h2 ∘ h1

# (4x+3) ∘ (2x+1) = 8x+7
b3 = Point(7,7)
m3 = TensorValue(8,0,0,8)
h3 = AffineMap(m3,b3)

h1inv = inverse_map(h1)
h2inv = inverse_map(h2)
hinv = inverse_map(h)
hinvinv = inverse_map(hinv)

function fields_are_equal(f::Field, g::Field)
    local xs = Point{2,Float64}[Point(0,0), Point(1,0), Point(0,1)]
    # all(f(x) ≈ g(x) for x in xs)
    all(f(xs) ≈ g(xs))
end

@test fields_are_equal(id ∘ id, id)

@test fields_are_equal(h, h3)

@test fields_are_equal(id ∘ h, h)
@test fields_are_equal(h ∘ id, h)

@test fields_are_equal(hinv ∘ h, id)
@test fields_are_equal(h ∘ hinv, id)

@test fields_are_equal(hinvinv, h)

@test fields_are_equal(hinv, h1inv ∘ h2inv)

end # module
