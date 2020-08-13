module CellFieldsTests

using Test
using Gridap.Arrays
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

isa(cf1,Geometry.CellField)
@test RefStyle(cf1) == Val{true}()

indices = [3,2,3]
_cf1 = reindex(cf1,indices)
@test length(_cf1) == length(indices)
test_cell_field(_cf1,reindex(q,indices),reindex(r,indices),grad=reindex(g,indices))
isa(_cf1,Geometry.CellField)
@test RefStyle(_cf1) == Val{true}()

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

div_cfu = ∇⋅cfu
collect(evaluate(div_cfu,q))

curl_cfu = cross(∇,cfu)
collect(evaluate(curl_cfu,q))

epsi_cfu = ε(cfu)
collect(evaluate(epsi_cfu,q))

btrian = BoundaryTriangulation(model)

s = CompressedArray([Point{1,Float64}[(0.25,),(0.75,)]],get_cell_type(btrian))

bcf1 = restrict(cf1,btrian)
@test isa(bcf1,CellField)

nvec = get_normal_vector(btrian)
@test isa(nvec,CellField)

z = 2*bcf1 + nvec
@test isa(z,CellField)

flux1 = nvec⋅∇u
collect(evaluate(flux1,s))

flux2 = ∇u⋅nvec
collect(evaluate(flux2,s))

strian = SkeletonTriangulation(model)

scf1 = restrict(cf1,strian)

@test scf1.left === scf1.inward
@test scf1.right === scf1.outward

@test isa(scf1,SkeletonCellField)

j = jump(scf1)

s = CompressedArray([Point{1,Float64}[(0.25,),(0.75,)]],get_cell_type(strian))
r = fill([0.0,0.0],length(j))
test_cell_field(j,s,r)

m = mean(scf1)
r = [[0.5, 1.5], [2.0, 2.0], [2.5, 3.5], [2.0, 2.0]]
test_cell_field(m,s,r)

j = jump(∇(scf1))
r = fill(VectorValue{2,Float64}[(0,0),(0,0)],length(j))
test_cell_field(j,s,r)

scfu = restrict(cfu,strian)

grad_cfu = jump(outer(∇,scfu))
collect(evaluate(grad_cfu,s))

t_grad_cfu = jump(outer(scfu,∇))
collect(evaluate(t_grad_cfu,s))

div_cfu = mean(∇⋅scfu)
collect(evaluate(div_cfu,s))

curl_cfu = jump(cross(∇,scfu))
collect(evaluate(curl_cfu,s))

epsi_cfu = mean(ε(scfu))
collect(evaluate(epsi_cfu,s))

#using Gridap.Visualization
#
#writevtk(btrian,"btrian",cellfields=["bcf1"=>bcf1,"nvec"=>nvec,"z"=>z])

#using Gridap.Visualization
#
#writevtk(strian,"strian",cellfields=["j"=>j])

end # module
