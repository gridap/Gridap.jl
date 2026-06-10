module DiffOperatorsTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Fields: MockField, MockFieldArray
using Gridap.Fields: skew_symmetric_gradient, ShiftedNabla

using LinearAlgebra
using FillArrays
using StaticArrays: @MMatrix

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

s = ∇ - v
@test evaluate(outer(f,s),x) == transpose.(evaluate(outer(s,f),x))

l = 10
f = Fill(f,l)

@test Broadcasting(divergence)(f) == Broadcasting(Operation(tr))(Broadcasting(∇)(f))
@test Broadcasting(curl)(f) == Broadcasting(Operation(grad2curl))(Broadcasting(∇)(f))
@test Broadcasting(ε)(f) == Broadcasting(Operation(symmetric_part))(Broadcasting(∇)(f))

let p=p
  g_plus_pf = (g,f)->g+p⊗f
  @test evaluate(Broadcasting(∇+p)(f),x) == evaluate( Broadcasting(Operation(g_plus_pf))(Broadcasting(∇)(f),f)  ,x)
end

# Test automatic differentiation

u_scal(x) = x[1]^2 + x[2]
∇u_scal(x) = VectorValue( 2*x[1], one(x[2]) )
Δu_scal(x) = 2

u_vec(x) = VectorValue( x[1]^2 + x[2], 4*x[1] - x[2]^2 )
∇u_vec(x) = TensorValue( 2*x[1], one(x[2]), 4*one(x[1]), - 2*x[2] )
Δu_vec(x) = VectorValue( 2, -2 )
εu_vec(x) = SymTensorValue( 2*x[1], 0.5*(one(x[2])+4*one(x[1])), - 2*x[2] )

u_ten(x) = TensorValue( x[1]^2 + x[2], 4*x[1] - x[2]^2, -x[1]^2 - x[2], -4*x[1] + x[2]^2 )
∇u_ten(x) = ThirdOrderTensorValue( 2*x[1], one(x[2]), 4*one(x[1]), - 2*x[2],
                                  -2*x[1],-one(x[2]),-4*one(x[1]), + 2*x[2] )
Δu_ten(x) = TensorValue( 2, -2, -2, 2 )

u_ten23(x) = TensorValue{2,3}( x[1]^2 + x[2], 4*x[1] - x[2]^2,
                              -x[1]^2 - x[2],-4*x[1] + x[2]^2,
                               x[1]^2 + x[2], 4*x[1] - x[2]^2)
∇u_ten23(x) = ThirdOrderTensorValue{2,2,3}(
                                   2*x[1], one(x[2]), 4*one(x[1]), - 2*x[2],
                                  -2*x[1],-one(x[2]),-4*one(x[1]), + 2*x[2],
                                   2*x[1], one(x[2]), 4*one(x[1]), - 2*x[2],)
Δu_ten23(x) = TensorValue{2,3}( 2, -2, -2, 2, 2, -2)

u_sten(x)  = SymTensorValue( x[1]^2 + x[2], 4*x[1] - x[2]^2, -4*x[1] + x[2]^2 )
∇u_sten(x) = ThirdOrderTensorValue( 2*x[1],      one(x[2]), 4*one(x[1]), - 2*x[2],
                                    4*one(x[1]),-2*x[2],   -4*one(x[1]), + 2*x[2] )
Δu_sten(x) = SymTensorValue( 2, -2, 2 )

u_qten(x)  = SymTracelessTensorValue( x[1]^3 + 2x[2]^3, 5*x[1]^3 - 7x[2]^3)
∇u_qten(x) = ThirdOrderTensorValue( 3x[1]^2,  6x[2]^2, 15x[1]^2, -21x[2]^2,
                                   15x[1]^2,-21x[2]^2, -3x[1]^2,  -6x[2]^2)
Δu_qten(x) = SymTracelessTensorValue( 6x[1] + 12x[2], 30x[1] - 42x[2] )

u_ten3(x)  = ThirdOrderTensorValue{2,1,2}( x[1]^2 + x[2], 4*x[1] - x[2]^2, -x[1]^2 - x[2], -4*x[1] + x[2]^2 )
Δu_ten3(x) = ThirdOrderTensorValue{2,1,2}( 2, -2, -2, 2 )

xs = [ Point(1.,2.), Point(1.,1.), Point(2.,0.), Point(0.,3.), Point(-1.,3.)]
for x in xs
  @test  ∇(u_scal)(x) == ∇u_scal(x)
  @test  Δ(u_scal)(x) == Δu_scal(x)

  @test  ∇(u_vec)(x)  == ∇u_vec(x)
  @test  Δ(u_vec)(x)  == Δu_vec(x)
  @test  ε(u_vec)(x)  == εu_vec(x)

  @test  ∇(u_ten)(x)  == ∇u_ten(x)
  @test  Δ(u_ten)(x)  == Δu_ten(x)
  @test (∇⋅u_ten)(x)  == tr(∇u_ten(x))
  @test  Δ(u_ten)(x)  == (∇⋅∇u_ten)(x)

  @test  ∇(u_ten23)(x)  == ∇u_ten23(x)
  @test  Δ(u_ten23)(x)  == Δu_ten23(x)
  @test (∇⋅u_ten23)(x)  == tr(∇u_ten23(x))
  @test  Δ(u_ten23)(x)  == (∇⋅∇u_ten23)(x)

  @test  ∇(u_sten)(x) == ∇u_sten(x)
  @test  Δ(u_sten)(x) == Δu_sten(x)
  @test (∇⋅u_sten)(x) == tr(∇u_sten(x))
 #@test  Δ(u_sten)(x) == (∇⋅∇u_sten)(x) # Would work if ∇(u_sten) is implemented to return a SymThirdOrderTensorValue
  @test  get_array(Δ(u_sten)(x)) == get_array((∇⋅∇u_sten)(x))

  @test  ∇(u_qten)(x) == ∇u_qten(x)
  @test  Δ(u_qten)(x) == Δu_qten(x)
  @test (∇⋅u_qten)(x) == tr(∇u_qten(x))
 #@test  Δ(u_qten)(x) == (∇⋅∇u_qten)(x)
  @test  get_array(Δ(u_qten)(x)) == get_array((∇⋅∇u_qten)(x))

  @test  Δ(u_ten3)(x)  == Δu_ten3(x)
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

let x0=VectorValue(-0.3, 0.5), A=TensorValue(0.0, -1.0, 1.0, 0.0)
  u(x) = (x - x0) ⋅ A
  ∇u(x) = A
  Δu(x) = VectorValue(0.0, 0.0)

  for x in xs
    @test ∇(u)(x) == ∇u(x)
    @test Δ(u)(x) == Δu(x)
    @test Δ(u)(x) == (∇⋅∇u)(x)
  end
end

let x0=VectorValue(1.0, 0.5)
  u(x) = exp(-(x - x0) ⋅ (x - x0))
  ∇u(x) = -2(x - x0) * u(x)
  Δu(x) = 4u(x) * ((x - x0) ⋅ (x - x0) - 1)

  for x in xs
    @test ∇(u)(x) == ∇u(x)
    @test Δ(u)(x) == Δu(x)
    @test Δ(u)(x) == (∇⋅∇u)(x)
  end
end

u(x) = VectorValue( x[1]^2 + 2*x[3]^2, -x[1]^2, -x[2]^2 + x[3]^2 )
∇u(x) = TensorValue( 2*x[1],0,4*x[3],  -2*x[1],0,0,  0,-2*x[2],2*x[3] )
Δu(x) = VectorValue( 6, -2, 0 )
εu(x) = symmetric_part(∇u(x))
νu(x) = TensorValues.skew_symmetric_part(∇u(x))

xs = [ Point(1.,1.,2.0), Point(2.,0.,1.), Point(0.,3.,0.), Point(-1.,3.,2.)]
for x in xs
  @test ∇(u)(x) == ∇u(x)
  @test (∇⋅u)(x) == tr(∇u(x))
  @test (∇×u)(x) == grad2curl(∇u(x))
  @test Δ(u)(x) == Δu(x)
  @test ε(u)(x) == εu(x)
  @test Fields.skew_symmetric_gradient(u)(x) == νu(x)
  @test Δ(u)(x) == (∇⋅∇u)(x)
  @test (∇⋅εu)(x) == VectorValue(4.,-1.,1.)
end

let v=VectorValue(0.5, -1.0, 2.0), x0=VectorValue(0.3, 1.0, -0.2)
  u(x) = (x - x0) × v + v ⋅ (x ⊗ x)
  function ∇u(x)
      a = TensorValue(@MMatrix [0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0] )
      b =  v ⊗ x + v ⋅ x * one(TensorValue{3,3,Float64})
      a+b
  end
  Δu(x) = 2v
  εu(x) = symmetric_part(∇u(x))
  νu(x) = TensorValues.skew_symmetric_part(∇u(x))

  xs = [ Point(1.,1.,2.0), Point(2.,0.,1.), Point(0.,3.,0.), Point(-1.,3.,2.)]
  x = first(xs)
  for x in xs
    @test ∇(u)(x) == ∇u(x)
    @test (∇⋅u)(x) == tr(∇u(x))
    @test (∇×u)(x) == grad2curl(∇u(x))
    @test Δ(u)(x) == Δu(x)
    @test ε(u)(x) == εu(x)
    @test Fields.skew_symmetric_gradient(u)(x) == νu(x)
    @test Δ(u)(x) == (∇⋅∇u)(x)
  end
end

# FieldArray diff ops

same_type(x,y) = typeof(x) == typeof(y)

x = [ Point(0.,0.), Point(1.,0.), Point(0.,1.), Point(1.,1.) ]
f = MockFieldArray( [ VectorValue(1.0,2.0), VectorValue(3.0,4.0), VectorValue(5.0,6.0) ] )

@test same_type(evaluate(Broadcasting(gradient),f), return_value(Broadcasting(gradient),f))
@test same_type(evaluate(Broadcasting(divergence),f), return_value(Broadcasting(divergence),f))
@test same_type(evaluate(Broadcasting(curl),f), return_value(Broadcasting(curl),f))
@test same_type(evaluate(Broadcasting(symmetric_gradient),f), return_value(Broadcasting(symmetric_gradient),f))
@test same_type(evaluate(Broadcasting(laplacian),f), return_value(Broadcasting(laplacian),f))
@test same_type(evaluate(Broadcasting(skew_symmetric_gradient),f), return_value(Broadcasting(skew_symmetric_gradient),f))

end # module
