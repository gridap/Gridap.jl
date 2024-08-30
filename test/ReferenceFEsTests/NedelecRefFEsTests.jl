module NedelecRefFEsTest

using Test
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Fields
using Gridap.TensorValues
using Gridap.Fields: MockField
using Gridap.ReferenceFEs

p = QUAD
D = num_dims(QUAD)
et = Float64
order = 0

reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test num_terms(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 0
@test num_dofs(reffe) == 4

@test Conformity(reffe) == CurlConformity()

p = QUAD
D = num_dims(QUAD)
et = Float64
order = 1

reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test num_terms(get_prebasis(reffe)) == 12
@test num_dofs(reffe) == 12
@test get_order(get_prebasis(reffe)) == 1

prebasis = get_prebasis(reffe)
dof_basis = get_dof_basis(reffe)

v = VectorValue(3.0,0.0)
field = GenericField(x->v*x[1])

cache = return_cache(dof_basis,field)
r = evaluate!(cache, dof_basis, field)
test_dof_array(dof_basis,field,r)

cache = return_cache(dof_basis,prebasis)
r = evaluate!(cache, dof_basis, prebasis)
test_dof_array(dof_basis,prebasis,r)

p = TET
D = num_dims(p)
et = Float64
order = 0

reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test num_terms(get_prebasis(reffe)) == 6
@test get_order(get_prebasis(reffe)) == 0
@test num_dofs(reffe) == 6
@test Conformity(reffe) == CurlConformity()

#using Gridap.Geometry
#p = TRI
#D = num_dims(p)
#et = Float64
#order = 0
#reffe = NedelecRefFE(et,p,order)
#grid = compute_reference_grid(p,4)
#x = get_node_coordinates(grid)
#shapes = get_prebasis(reffe)
#shapes = Broadcasting(∇)(shapes)
#ux = evaluate(shapes,x)
#for i in 8#1:size(ux,1)
#  @show i
#  @show x[i]
#  display(ux[i,:])
#end

#using Gridap.Geometry
#using Gridap.Visualization
#grid = compute_reference_grid(p,10)
#x = get_node_coordinates(grid)
#shapes = get_shapefuns(reffe)
#gshapes = Broadcasting(∇)(shapes)
#ux = evaluate(shapes,x)
#gux = evaluate(gshapes,x)
#ndat = ["s$i"=>ux[:,i] for i in 1:num_dofs(reffe)]
#gndat = ["g$i"=>gux[:,i] for i in 1:num_dofs(reffe)]
#writevtk(grid,"nede_tet",nodaldata=vcat(ndat,gndat))

p = TRI
D = num_dims(p)
et = Float64
order = 0

reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test num_terms(get_prebasis(reffe)) == 3
@test get_order(get_prebasis(reffe)) == 0
@test num_dofs(reffe) == 3
@test Conformity(reffe) == CurlConformity()
dof_basis = get_dof_basis(reffe)

#using Gridap.Geometry
#using Gridap.Visualization
#grid = compute_reference_grid(p,10)
#x = get_node_coordinates(grid)
#shapes = get_shapefuns(reffe)
#gshapes = Broadcasting(∇)(shapes)
#ux = evaluate(shapes,x)
#gux = evaluate(gshapes,x)
#ndat = ["s$i"=>ux[:,i] for i in 1:num_dofs(reffe)]
#gndat = ["g$i"=>gux[:,i] for i in 1:num_dofs(reffe)]
#writevtk(grid,"nede_tri",nodaldata=vcat(ndat,gndat))

# Factory function
reffe = ReferenceFE(QUAD,nedelec,0)
@test num_terms(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 0
@test num_dofs(reffe) == 4
@test Conformity(reffe) == CurlConformity()

reffe = ReferenceFE(QUAD,nedelec,Float64,0)
@test num_terms(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 0
@test num_dofs(reffe) == 4
@test Conformity(reffe) == CurlConformity()

@test Conformity(reffe,:L2) == L2Conformity()
@test Conformity(reffe,:Hcurl) == CurlConformity()
@test Conformity(reffe,:HCurl) == CurlConformity()

@test Nedelec() == nedelec

p = TET
D = num_dims(p)
et = Float64
order = 1
reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test num_terms(get_prebasis(reffe)) == 20
@test get_order(get_prebasis(reffe)) == 1
@test num_dofs(reffe) == 20
@test Conformity(reffe) == CurlConformity()
dof_basis = get_dof_basis(reffe)

face_odofs_L2 = get_face_own_dofs(reffe,L2Conformity())

@test face_odofs_L2 == [Int64[], Int64[], Int64[], Int64[], 
                     Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], 
                    collect(1:20)]

face_odofs    = get_face_own_dofs(reffe)
face_cdofs    = get_face_dofs(reffe)

@test face_odofs == [Int64[], Int64[], Int64[], Int64[], 
                    [1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [11, 12], [13, 14], [15, 16], [17, 18], [19, 20], 
                    Int64[]]

@test face_cdofs == [Int64[], Int64[], Int64[], Int64[], 
                     [1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [11, 12], 
                     [1, 2, 3, 4, 5, 6, 13, 14], [1, 2, 7, 8, 9, 10, 15, 16], [3, 4, 7, 8, 11, 12, 17, 18], [5, 6, 9, 10, 11, 12, 19, 20], 
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]]

                     
#display(face_odofs)

using Gridap.Geometry
using Gridap.Visualization
grid = compute_reference_grid(p,10)
x = get_node_coordinates(grid)
shapes = get_shapefuns(reffe)
gshapes = Broadcasting(∇)(shapes)
ux = evaluate(shapes,x)
gux = evaluate(gshapes,x)
ndat = ["s$i"=>ux[:,i] for i in 1:num_dofs(reffe)]
gndat = ["g$i"=>gux[:,i] for i in 1:num_dofs(reffe)]
d = mktempdir()
f = joinpath(d, "nede_tet_1")
writevtk(grid,f,nodaldata=vcat(ndat,gndat))
rm(d,recursive=true)

end # module
