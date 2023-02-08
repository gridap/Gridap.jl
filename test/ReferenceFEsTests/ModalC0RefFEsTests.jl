module ModalC0RefFEsTests

using Test
using Gridap
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.Fields
using Gridap.CellData

# using BenchmarkTools

order = 1
p = QUAD
T = VectorValue{2,Float64}

m = ReferenceFE(p,modalC0,T,order)
l = ReferenceFE(p,lagrangian,T,order)

test_reference_fe(m)

@test num_dofs(m) == num_dofs(l)
@test Conformity(m) === Conformity(l)
@test get_face_own_dofs(m,Conformity(m)) == get_face_own_dofs(l,Conformity(l))
@test get_face_own_dofs_permutations(m,Conformity(m)) == get_face_own_dofs_permutations(l,Conformity(l))

# domain = (0,1,0,1)
# partition = (2,2)
# model = CartesianDiscreteModel(domain,partition)

function _test_function_interpolation(order,u,V)
  test_single_field_fe_space(V)
  uh = interpolate(u,V)
  Ω = Triangulation(model)
  degree = 2*order
  dΩ = Measure(Ω,degree)
  l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
  e = u - uh
  el2 = l2(e)
  @test el2 < 1.0e-9
  # writevtk(Ω,"results",nsubcells=20,cellfields=["uh"=>uh])
  # free_vals = zeros(num_free_dofs(V)); free_vals[7] = 1
  # uh = FEFunction(V,free_vals)
  # writevtk(Ω,"shape_7",nsubcells=20,cellfields=["s7"=>uh])
  # free_vals = zeros(num_free_dofs(V)); free_vals[9] = 1
  # uh = FEFunction(V,free_vals)
  # writevtk(Ω,"shape_9",nsubcells=20,cellfields=["s9"=>uh])
  # free_vals = zeros(num_free_dofs(V)); free_vals[11] = 1
  # uh = FEFunction(V,free_vals)
  # writevtk(Ω,"shape_11",nsubcells=20,cellfields=["s11"=>uh])
  # free_vals = zeros(num_free_dofs(V)); free_vals[13] = 1
  # uh = FEFunction(V,free_vals)
  # writevtk(Ω,"shape_13",nsubcells=20,cellfields=["s13"=>uh])
end

function test_function_interpolation(::Type{T},order,C,u) where T
  reffe = ReferenceFE(modalC0,T,order)
  V = FESpace(model,reffe,conformity=C)
  _test_function_interpolation(order,u,V)
end

function test_function_interpolation(::Type{T},order,C,u,bboxes,space) where T
  reffe = ReferenceFE(modalC0,T,order,bboxes,space=space)
  V = FESpace(model,reffe,conformity=C)
  _test_function_interpolation(order,u,V)
end

domain = (0,1)
partition = (1,)
model = CartesianDiscreteModel(domain,partition)
trian = Triangulation(model)

T = Float64; order = 3; C = :H1; u(x) = x[1]^3
bboxes = [Point(0.0),Point(1.0),Point(0.0),Point(1.0),Point(-3.0),Point(1.0)]
bboxes = CellPoint(fill(bboxes,1),trian,PhysicalDomain())
test_function_interpolation(T,order,C,u,bboxes.cell_ref_point,:Q)

domain = (0,1)
partition = (4,)
model = CartesianDiscreteModel(domain,partition)

T = Float64; order = 3; C = :H1; u(x) = x[1]^3
bboxes = [ [Point(0.0),Point(1.0),Point(0.0),Point(1.0),Point(-1.0),Point(3.0)],
           [Point(0.0),Point(1.0),Point(0.0),Point(1.0),Point(-1.5),Point(3.2)],
           [Point(0.0),Point(1.0),Point(0.0),Point(1.0),Point(-0.2),Point(1.0)],
           [Point(0.0),Point(1.0),Point(0.0),Point(1.0),Point(-1.2),Point(2.0)] ]
test_function_interpolation(T,order,C,u,bboxes,:Q)

domain = (0,4)
partition = (4,)
model = CartesianDiscreteModel(domain,partition)

T = Float64; order = 3; C = :H1; u(x) = x[1]^3
bboxes = [ [Point(0.0),Point(1.0),Point(0.0),Point(1.0),Point(-1.0),Point(3.0)],
           [Point(0.0),Point(1.0),Point(0.0),Point(1.0),Point(-1.5),Point(3.2)],
           [Point(0.0),Point(1.0),Point(0.0),Point(1.0),Point(-0.2),Point(1.0)],
           [Point(0.0),Point(1.0),Point(0.0),Point(1.0),Point(-1.2),Point(2.0)] ]
test_function_interpolation(T,order,C,u,bboxes,:Q)

domain = (0,1,0,1)
partition = (2,2,)
model = CartesianDiscreteModel(domain,partition)
trian = Triangulation(model)

T = Float64; order = 3; C = :H1; u(x) = (x[1]+x[2])^3
bboxes = reshape( [Point(0.0,0.0),Point(1.0,1.0),Point(0.0,0.0),
                   Point(1.0,1.0),Point(0.0,0.0),Point(1.0,1.0),
                   Point(0.0,0.0),Point(1.0,1.0),Point(-0.5,2.5),
                   Point(2.5,2.5),Point(-0.5,2.5),Point(2.5,2.5),
                   Point(-1.0,-1.0),Point(-1.0,1.25),Point(-1.0,-1.0),
                   Point(-1.0,1.25),Point(0.0,0.0),Point(1.0,1.0)], 18 )
bboxes = CellPoint(fill(bboxes,4),trian,PhysicalDomain())
test_function_interpolation(T,order,C,u,bboxes.cell_ref_point,:Q)

T = Float64; order = 5; C = :H1; u(x) = (x[1]+x[2])^5
bboxes = reshape( [Point(0.0,0.0),Point(1.0,1.0),Point(0.0,0.0),
                   Point(1.0,1.0),Point(0.0,0.0),Point(1.0,1.0),
                   Point(0.0,0.0),Point(1.0,1.0),Point(-0.5,2.5),
                   Point(2.5,2.5),Point(-0.5,2.5),Point(2.5,2.5),
                   Point(-1.0,-1.0),Point(-1.0,1.25),Point(-1.0,-1.0),
                   Point(-1.0,1.25),Point(0.0,0.0),Point(1.0,1.0)], 18 )
bboxes = CellPoint(fill(bboxes,4),trian,PhysicalDomain())
test_function_interpolation(T,order,C,u,bboxes.cell_ref_point,:S)

order = 1; T = Float64; C = :L2; u(x) = x[1]+x[2]
test_function_interpolation(T,order,C,u)

order = 1; T = VectorValue{2,Float64}; C = :H1
u(x) = VectorValue(x[1]+x[2],x[2])
test_function_interpolation(T,order,C,u)

domain = (0,1,0,1,0,1)
partition = (2,2,2)
model = CartesianDiscreteModel(domain,partition)

order = 1; T = Float64; C = :H1; u(x) = x[1]+x[2]+x[3]
test_function_interpolation(T,order,C,u)

# Inspect operator matrix to check if L2-scalar product of
# gradients of bubble functions satisfy Kronecker's delta
# domain = (0,1)
# partition = (1)
# model = CartesianDiscreteModel(domain,partition)
# order = 6; T = Float64; C = :H1;
# reffe = ReferenceFE(modalC0,T,order)
# V = FESpace(model,reffe,conformity=C)
# Ω = Triangulation(model)
# degree = 2*order
# dΩ = LebesgueMeasure(Ω,degree)
# a(u,v) = ∫( ∇(v)⊙∇(u) )*dΩ
# b(v) = 0.0
# op = AffineFEOperator(bboxes,V,V)

end # module
