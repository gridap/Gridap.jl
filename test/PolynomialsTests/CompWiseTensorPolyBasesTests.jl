module CompWiseTensorPolyBasisTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials

PT = Bernstein
T = Float64

D = 2
x = [Point(1.,0.), Point(.0,.5), Point(.2,.3)]

orders = Int[ 3 0; 1 0]
V = VectorValue{2,T}
b = CompWiseTensorPolyBasis{D}(PT,V, orders)

@test Tuple.(b.comp_terms[1]) == [ (1, 1), (2, 1), (3, 1), (4, 1)]
@test Tuple.(b.comp_terms[2]) == [ (1, 1), (2, 1), ]

x₁ = [Point(1.), Point(.0), Point(.2)]
b_1x = BernsteinBasis(Val(1),T,3)
b_2x = BernsteinBasis(Val(1),T,1)
Gb_1x = Broadcasting(∇)(b_1x)
Gb_2x = Broadcasting(∇)(b_2x)
Hb_1x = Broadcasting(∇∇)(b_1x)
Hb_2x = Broadcasting(∇∇)(b_2x)

bx  = hcat( [ bi .* VectorValue(1.,0) for bi in evaluate(b_1x, x₁) ],
            [ bi .* VectorValue(0.,1) for bi in evaluate(b_2x, x₁) ])
Gbx = hcat( [ getindex.(bi,1) .* TensorValue(1.,0,0,0) for bi in evaluate(Gb_1x, x₁) ],
            [ getindex.(bi,1) .* TensorValue(0.,0,1,0) for bi in evaluate(Gb_2x, x₁) ])
Hbx = hcat( [ getindex.(bi,1) .* ThirdOrderTensorValue(1.,0,0,0,0,0,0,0) for bi in evaluate(Hb_1x, x₁) ],
            [ getindex.(bi,1) .* ThirdOrderTensorValue(0.,0,0,0,1,0,0,0) for bi in evaluate(Hb_2x, x₁) ])

test_field_array(b,x,bx,≈, grad=Gbx, gradgrad=Hbx)

# Dependent components

V = SymTracelessTensorValue{2,T}
b = CompWiseTensorPolyBasis{D}(PT,V, orders)
evaluate(b,x)
evaluate(Broadcasting(∇)(b),x)

@test Tuple.(b.comp_terms[1]) == [ (1, 1), (2, 1), (3, 1), (4, 1)]
@test Tuple.(b.comp_terms[2]) == [ (1, 1), (2, 1), ]

x₁ = [Point(1.), Point(.0), Point(.2)]
b_1x = BernsteinBasis(Val(1),T,3)
b_2x = BernsteinBasis(Val(1),T,1)
Gb_1x = Broadcasting(∇)(b_1x)
Gb_2x = Broadcasting(∇)(b_2x)

bx  = hcat( [ bi .* V(1.,0) for bi in evaluate(b_1x, x₁) ],
            [ bi .* V(0.,1) for bi in evaluate(b_2x, x₁) ])
Gbx = hcat( [ getindex.(bi,1) .* VectorValue(1.,0)⊗V(1.,0) for bi in evaluate(Gb_1x, x₁) ],
            [ getindex.(bi,1) .* VectorValue(1.,0)⊗V(0.,1) for bi in evaluate(Gb_2x, x₁) ])

test_field_array(b,x,bx,≈, grad=Gbx)


end # module
