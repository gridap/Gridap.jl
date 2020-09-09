module CellFieldsTests

using Test
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.TensorValues
using FillArrays
using LinearAlgebra
using Gridap.CellData

import Gridap.Fields: ∇

domain = (0,4,0,4)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

trian = get_triangulation(model)
q2x = get_cell_map(trian)

qi = [Point(0.0,0.0),]
l = num_cells(trian)
q = Fill(qi,l)

fun(x) = x[1]
∇fun(x) = VectorValue(1,0)
∇(::typeof(fun)) = ∇fun

cf1 = CellField(fun,trian)
r = [[0.0], [2.0], [0.0], [2.0]]
g = Vector{VectorValue{2,Float64}}[[(1, 0)], [(1, 0)], [(1, 0)], [(1, 0)]]
test_cell_field(cf1∘q2x,q,r,grad=g)

isa(cf1,Geometry.CellField)

indices = [3,2,3]
_cf1 = reindex(cf1∘q2x,indices)
@test length(_cf1) == length(indices)
test_cell_field(_cf1,reindex(q,indices),reindex(r,indices),grad=reindex(g,indices))
isa(_cf1,Geometry.CellField)

cf2 = cf1 + 4
r2 = map( (i) -> i.+4, r)
test_cell_field(cf2∘q2x,q,r2)

cf3 = 2*cf2
r3 = map( (i) -> 2*i, r2)
test_cell_field(cf3∘q2x,q,r3)

cf4 = -cf1
r4 = -r
g4 = -g
test_cell_field(cf4∘q2x,q,r4,grad=g4)

u(x) = VectorValue(x[1],x[2])
∇u(x) = TensorValue(1,0,0,1)
∇(::typeof(u)) = ∇u

cfu = CellField(u,trian)
collect(evaluate(cfu∘q2x,q))

grad_cfu = outer(∇,cfu)
collect(evaluate(grad_cfu∘q2x,q))

t_grad_cfu = outer(cfu,∇)
collect(evaluate(t_grad_cfu∘q2x,q))

div_cfu = ∇⋅cfu
collect(evaluate(div_cfu∘q2x,q))

curl_cfu = cross(∇,cfu)
collect(evaluate(curl_cfu∘q2x,q))

epsi_cfu = ε(cfu)
collect(evaluate(epsi_cfu∘q2x,q))

btrian = BoundaryTriangulation(model)

s = CompressedArray([Point{1,Float64}[(0.25,),(0.75,)]],get_cell_type(btrian))

s2x = get_cell_map(btrian)

bcf1 = cf1∘s2x
@test isa(bcf1,CellField)

nvec = get_normal_vector(btrian)
@test isa(nvec,CellField)

z = 2*cf1 + nvec
@test isa(z,CellField)

flux1 = nvec⋅∇u
collect(evaluate(flux1∘s2x,s))

flux2 = ∇u⋅nvec
collect(evaluate(flux2∘s2x,s))

strian = SkeletonTriangulation(model)

z2x = get_cell_map(strian)

j = jump(cf1)

s = CompressedArray([Point{1,Float64}[(0.25,),(0.75,)]],get_cell_type(strian))
r = fill([0.0,0.0],length(j))
test_cell_field(j∘z2x,s,r)

m = mean(cf1)
r = [[0.5, 1.5], [2.0, 2.0], [2.5, 3.5], [2.0, 2.0]]
test_cell_field(m∘z2x,s,r)

j = jump(∇(cf1))
r = fill(VectorValue{2,Float64}[(0,0),(0,0)],length(j))
test_cell_field(j∘z2x,s,r)

grad_cfu = jump(outer(∇,cfu))
collect(evaluate(grad_cfu∘z2x,s))

t_grad_cfu = jump(outer(cfu,∇))
collect(evaluate(t_grad_cfu∘z2x,s))

div_cfu = mean(∇⋅cfu)
collect(evaluate(div_cfu∘z2x,s))

curl_cfu = jump(cross(∇,cfu))
collect(evaluate(curl_cfu∘z2x,s))

epsi_cfu = mean(ε(cfu))
collect(evaluate(epsi_cfu∘z2x,s))

#using Gridap.Visualization
#
#writevtk(btrian,"btrian",cellfields=["bcf1"=>bcf1,"nvec"=>nvec,"z"=>z])

#using Gridap.Visualization
#
#writevtk(strian,"strian",cellfields=["j"=>j])

end # module
