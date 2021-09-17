module SkeletonTriangulationsTests

using Test
using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: IN, OUT
using Gridap.Geometry: DiscreteModelMock

model = DiscreteModelMock()

strian = SkeletonTriangulation(model)
test_triangulation(strian)
@test get_background_model(strian) === model
glue = get_glue(strian,Val(1))
@test glue.tface_to_mface === strian.plus.glue.face_to_bgface
glue = get_glue(strian,Val(2))
@test glue.plus.tface_to_mface === strian.plus.glue.face_to_cell
@test glue.minus.tface_to_mface === strian.minus.glue.face_to_cell

domain = (0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)

strian = SkeletonTriangulation(model)
test_triangulation(strian)
@test get_background_model(strian) === model
glue = get_glue(strian,Val(2))
@test glue.tface_to_mface === strian.plus.glue.face_to_bgface
glue = get_glue(strian,Val(3))
@test glue.plus.tface_to_mface === strian.plus.glue.face_to_cell
@test glue.minus.tface_to_mface === strian.minus.glue.face_to_cell

model = DiscreteModelMock()

trian = Triangulation(model)

cell_to_is_left = [true,true,false,true,false]

itrian = InterfaceTriangulation(model,cell_to_is_left)

ltrian = itrian.plus
rtrian = itrian.minus
@test itrian.⁺ === ltrian
@test itrian.⁻ === rtrian

ni = get_facet_normal(itrian)
nl = get_facet_normal(ltrian)
nr = get_facet_normal(rtrian)
@test isa(ni,SkeletonPair)

#using Gridap.Visualization
#
#writevtk(trian,"trian")
#writevtk(itrian,"itrian",nsubcells=10,cellfields=["ni"=>ni,"nl"=>nl,"nr"=>nr])

domain = (0,1,0,1)
partition = (10,10)
model = CartesianDiscreteModel(domain,partition)

cell_to_inout = zeros(Int8,num_cells(model))
cell_to_inout[1:13] .= OUT
cell_to_inout[14:34] .= IN

itrian = InterfaceTriangulation(model,cell_to_inout)
@test num_cells(itrian) == 11

itrian = InterfaceTriangulation(model,1:13,14:34)
@test num_cells(itrian) == 11
Ω_in = Triangulation(model,findall(i->i==IN,cell_to_inout))
Ω_out = Triangulation(model,findall(i->i==OUT,cell_to_inout))
itrian = InterfaceTriangulation(Ω_in,Ω_out)

#ltrian = get_left_boundary(itrian)
#rtrian = get_right_boundary(itrian)
#ni = get_normal_vector(itrian)
#nl = get_normal_vector(ltrian)
#nr = get_normal_vector(rtrian)
#trian = Triangulation(model)
#
#using Gridap.Visualization
#writevtk(trian,"trian",celldata=["inout"=>cell_to_inout])
#writevtk(itrian,"itrian",nsubcells=10,cellfields=["ni"=>ni,"nl"=>nl,"nr"=>nr])

reffe = LagrangianRefFE(Float64,QUAD,(2,2))
conf = CDConformity((CONT,DISC))
face_own_dofs = get_face_own_dofs(reffe,conf)
strian = SkeletonTriangulation(model,reffe,face_own_dofs)
test_triangulation(strian)
ns = get_facet_normal(strian)
@test length(ns.⁺) == num_cells(strian)

#using Gridap.Visualization
#writevtk(strian,"strian",cellfields=["normal"=>ns])

domain = (0,1,0,1,0,1)
partition = (3,3,3)
oldmodel = CartesianDiscreteModel(domain,partition)
oldstrian = SkeletonTriangulation(oldmodel)

end # module
