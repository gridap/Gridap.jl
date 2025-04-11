module SkeletonTriangulationsTests

using Test
using Gridap
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

ids = [1,3,4]
vtrian = view(strian,ids)
test_triangulation(vtrian)
@test num_cells(vtrian) == 3
sglue = get_glue(strian,Val(2))
vglue = get_glue(vtrian,Val(2))
@test vglue.tface_to_mface == sglue.tface_to_mface[ids]
sglue = get_glue(strian,Val(3))
vglue = get_glue(vtrian,Val(3))
@test vglue.plus.tface_to_mface == sglue.plus.tface_to_mface[ids]
@test vglue.minus.tface_to_mface == sglue.minus.tface_to_mface[ids]
vn = get_facet_normal(vtrian)
@test isa(vn,SkeletonPair)
@test isa(vn.plus,AbstractArray)
@test isa(vn.minus,AbstractArray)

Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
Λ = SkeletonTriangulation(Γ)
normal = get_normal_vector(Λ)
@test Λ.rtrian === Γ
@test isa(Λ.dtrian,SkeletonTriangulation)
glue = get_glue(Λ,Val(3))
glue = get_glue(Λ,Val(2))
@test is_change_possible(Ω,Λ)
@test is_change_possible(Γ,Λ)

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

ti = get_tangent_vector(itrian)
tl = get_tangent_vector(ltrian)
tr = get_tangent_vector(rtrian)

@test ti isa SkeletonPair
@test tl isa Gridap.CellData.GenericCellField
@test tr isa Gridap.CellData.GenericCellField 

reffe = LagrangianRefFE(Float64,QUAD,(2,2))
conf = CDConformity((CONT,DISC))
face_own_dofs = get_face_own_dofs(reffe,conf)
strian = SkeletonTriangulation(model,reffe,face_own_dofs)
test_triangulation(strian)
ns = get_facet_normal(strian)
ts = get_edge_tangent(strian)
@test length(ns.⁺) == num_cells(strian)
@test length(ts.⁺) == num_cells(strian)
# Test orthogonality
should_be_zero = lazy_map(Broadcasting(Operation(dot)), ns.plus, ts.plus)
ndot_t_evaluated = evaluate(should_be_zero, VectorValue.(0:0.05:1))
@test maximum(ndot_t_evaluated) < 1e-15

#using Gridap.Visualization
#writevtk(strian,"strian",cellfields=["normal"=>ns])

domain = (0,1,0,1,0,1)
partition = (3,3,3)
oldmodel = CartesianDiscreteModel(domain,partition)
oldstrian = SkeletonTriangulation(oldmodel)

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
Ω1 = Interior(model,[2])
Ω2 = Interior(model,[4])
Γ12 = Interface(Ω1,Ω2)
@test num_cells(Γ12) == 1


end # module
