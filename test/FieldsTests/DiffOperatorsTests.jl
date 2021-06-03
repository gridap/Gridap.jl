module DiffOperatorsTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Fields: MockField
using LinearAlgebra
using FillArrays

np = 4
p = Point(1,2)
x = fill(p,np)

v = VectorValue(3.0,2.0)
f = MockField(v)

@test ∇(f) == gradient(f)

@test divergence(f) == tr(gradient(f))

@test curl(f) == grad2curl(gradient(f))

@test ∇⋅f == divergence(f)

@test cross(∇,f) == curl(f)

@test ∇×f == curl(f)

@test outer(∇,f) == ∇(f)

@test ∇⊗f == ∇(f)

@test outer(f,∇) == transpose(∇(f))

@test f⊗∇ == transpose(∇(f))

@test ε(f) == symmetric_part(gradient(f))

@test Δ(f) == ∇⋅∇(f)

@test (∇.*f)(x) != nothing

@test ((∇+p).*f)(x) != nothing

@test ((∇+p)(f))(x) == (∇(f) + p⊗f)(x)

g(x) = 2*x[2]

@test ((∇+p)(g))(x) == (∇(GenericField(g)) + p⊗GenericField(g))(x)

@test (∇+p)⋅f != nothing
@test (∇+p)×f != nothing
@test (∇+p)⊗f != nothing
@test f⊗(∇+p) != nothing

l = 10
f = Fill(f,l)

@test Broadcasting(divergence)(f) == Broadcasting(Operation(tr))(Broadcasting(∇)(f))
@test Broadcasting(curl)(f) == Broadcasting(Operation(grad2curl))(Broadcasting(∇)(f))
@test Broadcasting(ε)(f) == Broadcasting(Operation(symmetric_part))(Broadcasting(∇)(f))

@test evaluate(Broadcasting(∇+p)(f),x) == evaluate( Broadcasting(Operation((g,f)->g+p⊗f))(Broadcasting(∇)(f),f)  ,x)

# Test automatic differentiation

u_scal(x) = x[1]^2 + x[2]
∇u_scal(x) = VectorValue( 2*x[1], one(x[2]) )
Δu_scal(x) = 2

u_vec(x) = VectorValue( x[1]^2 + x[2], 4*x[1] - x[2]^2 )
∇u_vec(x) = TensorValue( 2*x[1], one(x[2]), 4*one(x[1]), - 2*x[2] )
Δu_vec(x) = VectorValue( 2, -2 )
εu_vec(x) = SymTensorValue( 2*x[1], 0.5*(one(x[2])+4*one(x[1])), - 2*x[2] )

xs = [ Point(1.,1.), Point(2.,0.), Point(0.,3.), Point(-1.,3.)]
for x in xs
  @test ∇(u_scal)(x) == ∇u_scal(x)
  @test Δ(u_scal)(x) == Δu_scal(x)
  @test ∇(u_vec)(x) == ∇u_vec(x)
  @test Δ(u_vec)(x) == Δu_vec(x)
  @test ε(u_vec)(x) == εu_vec(x)
end

u(x) = VectorValue( x[1]^2 + 2*x[2]^2, -x[1]^2 )
∇u(x) = TensorValue( 2*x[1], 4*x[2], -2*x[1], zero(x[1]) )
Δu(x) = VectorValue( 6, -2 )
εu(x) = SymTensorValue( 2*x[1], 2*x[2]-x[1], zero(x[1]) )

for x in xs
  @test (∇⋅u)(x) == tr(∇u(x)) 
  @test (∇×u)(x) == grad2curl(∇u(x))
  @test Δ(u)(x) == Δu(x)
  @test ε(u)(x) == εu(x)
  @test ∇(u)(x) == ∇u(x)
  @test Δ(u)(x) == (∇⋅∇u)(x)
  @test (∇⋅εu)(x) == VectorValue(4.,-1.)
end

u(x) = VectorValue( x[1]^2 + 2*x[2]^2, 0 )
∇u(x) = TensorValue( 2*x[1], 4*x[2], 0, 0 )
Δu(x) = VectorValue( 6, 0 )
εu(x) = SymTensorValue( 2*x[1], 2*x[2], 0 )

for x in xs
  @test (∇⋅u)(x) == tr(∇u(x)) 
  @test (∇×u)(x) == grad2curl(∇u(x))
  @test Δ(u)(x) == Δu(x)
  @test ε(u)(x) == εu(x)
  @test Δ(u)(x) == (∇⋅∇u)(x)
  @test (∇⋅εu)(x) == VectorValue(4.,0.)
end

u(x) = VectorValue( x[1]^2 + 2*x[3]^2, -x[1]^2, -x[2]^2 + x[3]^2 )
∇u(x) = TensorValue( 2*x[1],0,4*x[3],  -2*x[1],0,0,  0,-2*x[2],2*x[3] )
Δu(x) = VectorValue( 6, -2, 0 )
εu(x) = symmetric_part(∇u(x))

xs = [ Point(1.,1.,2.0), Point(2.,0.,1.), Point(0.,3.,0.), Point(-1.,3.,2.)]
for x in xs
  @test ∇(u)(x) == ∇u(x)
  @test (∇⋅u)(x) == tr(∇u(x)) 
  @test (∇×u)(x) == grad2curl(∇u(x))
  @test Δ(u)(x) == Δu(x)
  @test ε(u)(x) == εu(x)
  @test Δ(u)(x) == (∇⋅∇u)(x)
  @test (∇⋅εu)(x) == VectorValue(4.,-1.,1.)
end

end # module
