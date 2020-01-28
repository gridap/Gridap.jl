module DiffOperatorsTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Fields: MockField
using LinearAlgebra
using FillArrays

np = 4
p = Point(1,2)
x = fill(p,np)

v = 3.0
d = 2
_f = MockField{d}(v)

l = 10
_af = Fill(_f,l)

for f in (_f,_af)

  @test ∇(f) == gradient(f)
  
  @test divergence(f) == tr(gradient(f))
  
  @test curl(f) == grad2curl(gradient(f))
  
  @test ∇*f == divergence(f)
  
  @test cross(∇,f) == curl(f)
  
  @test outer(∇,f) == ∇(f)
  
  @test outer(f,∇) == transpose(∇(f))

  @test ε(f) == symmetic_part(gradient(f))

  @test Δ(f) == ∇*∇(f)

end

end # module
