module NedelecRefFEsTest

using Test
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Fields
using Gridap.TensorValues
using Gridap.Fields: MockField
using Gridap.ReferenceFEs


@test Nedelec{1}() == nedelec
@test Nedelec{1}() == nedelec1
@test Nedelec{2}() == nedelec2

p = QUAD
D = num_dims(QUAD)
et = Float64
order = 0

reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test length(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 1
@test num_dofs(reffe) == 4

@test Conformity(reffe) == CurlConformity()

@test_warn "falling back to `change_dof=false`" NedelecRefFE(et,p,order; change_dof=true, poly_type=Monomial)

p = QUAD
D = num_dims(QUAD)
et = Float64
order = 1

reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test length(get_prebasis(reffe)) == 12
@test num_dofs(reffe) == 12
@test get_order(get_prebasis(reffe)) == 2

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

order = 3
reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)

p = HEX
D = num_dims(HEX)
et = Float64
order = 2

reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test length(get_prebasis(reffe)) == 144
@test num_dofs(reffe) == 144
@test get_order(get_prebasis(reffe)) == 3

prebasis = get_prebasis(reffe)
dof_basis = get_dof_basis(reffe)

v = VectorValue(3.0,1.0,-2.)
field = GenericField(x->v*x[1]*x[2]*x[3])

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
@test length(get_prebasis(reffe)) == 6
@test get_order(get_prebasis(reffe)) == 1
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
@test length(get_prebasis(reffe)) == 3
@test get_order(get_prebasis(reffe)) == 1
@test num_dofs(reffe) == 3
@test Conformity(reffe) == CurlConformity()
dof_basis = get_dof_basis(reffe)

order = 3

reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test length(get_prebasis(reffe)) == 24
@test get_order(get_prebasis(reffe)) == 4
@test num_dofs(reffe) == 24
@test Conformity(reffe) == CurlConformity()
dof_basis = get_dof_basis(reffe)

# tests for Nedelec elements of the second kind
@test_throws AssertionError NedelecRefFE(Float64, QUAD, 1; kind=2)
@test_throws AssertionError NedelecRefFE(Float64, HEX, 2; kind=2)
@test_throws AssertionError NedelecRefFE(Float64, TRI, 0; kind=2)

p = TRI
et = Float64
order = 1

reffe2 = NedelecRefFE(et,p,order; kind=2)
test_reference_fe(reffe2)
@test length(get_prebasis(reffe2)) == 6
@test get_order(get_prebasis(reffe2)) == 1
@test num_dofs(reffe2) == 6
@test Conformity(reffe2) == CurlConformity()
dof_basis = get_dof_basis(reffe2)

p = TRI
et = Float64
order = 2

reffe2 = NedelecRefFE(et,p,order; kind=2)
test_reference_fe(reffe2)
@test length(get_prebasis(reffe2)) == 12
@test get_order(get_prebasis(reffe2)) == 2
@test num_dofs(reffe2) == 12
@test Conformity(reffe2) == CurlConformity()
dof_basis = get_dof_basis(reffe2)


p = TET
et = Float64
order = 3

reffe2 = NedelecRefFE(et,p,order; kind=2)
test_reference_fe(reffe2)
@test length(get_prebasis(reffe2)) == 60
@test get_order(get_prebasis(reffe2)) == 3
@test num_dofs(reffe2) == 60
@test Conformity(reffe2) == CurlConformity()
dof_basis = get_dof_basis(reffe2)


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
@test length(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 1
@test num_dofs(reffe) == 4
@test Conformity(reffe) == CurlConformity()

@test_warn "falling back to `change_dof=false`" ReferenceFE(QUAD,nedelec,0; poly_type=Monomial)

reffe = ReferenceFE(QUAD,nedelec,Float64,0)
@test length(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 1
@test num_dofs(reffe) == 4
@test Conformity(reffe) == CurlConformity()

@test Conformity(reffe,:L2) == L2Conformity()
@test Conformity(reffe,:Hcurl) == CurlConformity()
@test Conformity(reffe,:HCurl) == CurlConformity()


p = TET
D = num_dims(p)
et = Float64
order = 1
reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test length(get_prebasis(reffe)) == 20
@test get_order(get_prebasis(reffe)) == 2
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

p = TET
D = num_dims(p)
et = Float64
order = 3
reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test length(get_prebasis(reffe)) == 84
@test get_order(get_prebasis(reffe)) == 4
@test num_dofs(reffe) == 84
@test Conformity(reffe) == CurlConformity()
dof_basis = get_dof_basis(reffe)

face_odofs_L2 = get_face_own_dofs(reffe,L2Conformity())

@test face_odofs_L2 == [Int64[], Int64[], Int64[], Int64[],
                     Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[],
                    collect(1:84)]

face_odofs    = get_face_own_dofs(reffe)
face_cdofs    = get_face_dofs(reffe)

@test face_odofs == [Int64[], Int64[], Int64[], Int64[],
  [1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16], [17, 18, 19, 20], [21, 22, 23, 24],
  [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36],  # facet 123
  [37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48],  # facet 124
  [49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60],  # facet 134
  [61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72],  # facet 234
  [73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84]]  # cell

# Factory function: r=1, k=1
reffe = ReferenceFE(QUAD,nedelec,0)
@test reffe == ReferenceFE(QUAD,:Q⁻,1,1)

reffe = ReferenceFE(QUAD,nedelec,Float64,0)
@test reffe == ReferenceFE(QUAD,:Q⁻,1,1, Float64)

reffe = ReferenceFE(TRI,nedelec,0)
@test reffe == ReferenceFE(TRI,:P⁻,1,1)

reffe = ReferenceFE(HEX,nedelec,0)
@test reffe == ReferenceFE(HEX,:Q⁻,1,1)

reffe = ReferenceFE(TET,nedelec,0)
@test reffe == ReferenceFE(TET,:P⁻,1,1)

reffe2 = ReferenceFE(TRI,nedelec2,1)
@test reffe2 == ReferenceFE(TRI,:P,1,1)

reffe2 = ReferenceFE(TET,nedelec2,1)
@test reffe2 == ReferenceFE(TET,:P,1,1)

# Serendipity not implemented
@test_throws ErrorException ReferenceFE(QUAD,:S,1,1)
@test_throws ErrorException ReferenceFE(HEX, :S,1,1)

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

end # module
