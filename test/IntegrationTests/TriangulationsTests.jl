module TriangulationsTests

using Gridap
using Test
import Gridap: ∇

# CartesianTriangulation

partition = (2,4)
trian = CartesianTriangulation(partition)

test_triangulation(trian)

@test ncells(trian) == 8

coords = CellPoints(trian)

x4 = Point{2,Float64}[[0.0, -0.5], [1.0, -0.5], [0.0, 0.0], [1.0, 0.0]] 

@test coords[4] ≈ x4

# CellQuadrature

quad = CellQuadrature(trian,order=2)
q = coordinates(quad)
w = weights(quad)

# CellField

ufun(x::Point{2}) = 2*x[1]+x[2]

u = CellField(trian,ufun)
@test isa(u,CellField{2,Float64})

ufun_grad(x) = VectorValue(2.0,1.0)
∇(::typeof(ufun)) = ufun_grad

u = CellField(trian,ufun)
@test isa(u,CellField{2,Float64})

νfun(x::Point{2},u::Float64) = TensorValue(x[1], u*x[2], 0.0, u)

ν = CellField(trian,νfun,u)

@test isa(u,CellField{2,Float64})
@test isa(ν,CellField{2,TensorValue{2,Float64,4}})

#ν(u) = CellField(trian,νfun,u)
#w = ν(u) # julia nightly builds get stuck here

end # module
