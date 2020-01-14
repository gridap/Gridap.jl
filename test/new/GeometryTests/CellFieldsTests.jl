module CellFieldsTests

using Test
using Gridap.Geometry
using Gridap.Fields
using Gridap.TensorValues
using FillArrays
using LinearAlgebra

import Gridap.Fields: ∇

domain = (0,4,0,4)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

trian = get_triangulation(model)

qi = [Point(0.0,0.0),]
l = num_cells(trian)
q = Fill(qi,l)

fun(x) = x[1]
∇fun(x) = VectorValue(1,0)
∇(::typeof(fun)) = ∇fun

cf1 = CellField(fun,trian)
r = [[0.0], [2.0], [0.0], [2.0]]
g = Vector{VectorValue{2,Float64}}[[(1, 0)], [(1, 0)], [(1, 0)], [(1, 0)]]
test_cell_field(cf1,q,r,grad=g)

cf2 = cf1 + 4
r2 = map( (i) -> i.+4, r)
test_cell_field(cf2,q,r2)

cf3 = 2*cf2
r3 = map( (i) -> 2*i, r2)
test_cell_field(cf3,q,r3)

cf4 = -cf1
r4 = -r
g4 = -g
test_cell_field(cf4,q,r4,grad=g4)

u(x) = VectorValue(x[1],x[2])
∇u(x) = TensorValue(1,0,0,1)
∇(::typeof(u)) = ∇u

cfu = CellField(u,trian)
collect(evaluate(cfu,q))

grad_cfu = outer(∇,cfu)
collect(evaluate(grad_cfu,q))

t_grad_cfu = outer(cfu,∇)
collect(evaluate(t_grad_cfu,q))

div_cfu = ∇*cfu
collect(evaluate(div_cfu,q))

curl_cfu = cross(∇,cfu)
collect(evaluate(curl_cfu,q))

epsi_cfu = ε(cfu)
collect(evaluate(epsi_cfu,q))

btrian = BoundaryTriangulation(model)

bcf1 = restrict(cf1,btrian)
@test isa(bcf1,CellField)

nvec = get_normal_vector(btrian)
@test isa(nvec,CellField)

z = 2*bcf1 + nvec
@test isa(z,CellField)

#using Gridap.Visualization
#
#writevtk(btrian,"btrian",cellfields=["bcf1"=>bcf1,"nvec"=>nvec,"z"=>z])





end # module
