module AdaptedGeometryTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using FillArrays

for D = 1:3
  domain = Tuple(repeat([0,1],D))

  cart_model = CartesianDiscreteModel(domain,Tuple(fill(2,D)))
  model1 = refine(cart_model,4)
  model2 = refine(model1,Tuple(collect(1:D)))
  model3 = refine(model2,Tuple(fill(1,D)))
  test_discrete_model(model1)
  test_discrete_model(model2)
  test_discrete_model(model3)

  ctrian = Triangulation(cart_model)
  trian1 = Triangulation(model1)
  trian2 = Triangulation(model2)
  test_triangulation(trian1)
  test_triangulation(trian2)
  test_triangulation(trian1.trian)
  @test isa(trian1, AdaptedTriangulation)
  @test Gridap.Adaptivity.is_child(trian1,ctrian) == true
  @test Gridap.Adaptivity.is_child(ctrian,trian1) == false

  vtrian = view(trian1,[2,3,4])
  rtrian = Triangulation(trian1)

  # Get members
  fmodel = get_model(model1)
  cmodel = get_parent(model1)
  glue   = get_adaptivity_glue(model1)
  @test cmodel === cart_model
  @test fmodel === get_parent(model2)

  # Choosing targets
  ftrian = Triangulation(fmodel)
  ctrian = Triangulation(cmodel)
  @test best_target(trian1,ctrian) === trian1
  @test best_target(trian1,trian2) === trian2

  # Checking compatibility with other types of Triangulations
  t   = Triangulation(get_model(model1))
  rt  = Triangulation(model1)

  bt  = BoundaryTriangulation(model1)
  @test isa(bt,AdaptedTriangulation)
  test_triangulation(bt)
  bt_normal = get_normal_vector(bt)
  @test is_change_possible(t,bt)
  @test is_change_possible(rt,bt)
  @test !is_change_possible(bt,t)
  @test !is_change_possible(bt,rt)

  if D >= 2
    st  = SkeletonTriangulation(model1)
    @test isa(st,AdaptedTriangulation)
    test_triangulation(st)
    st_normal = get_normal_vector(st)
    @test is_change_possible(t,st)
    @test is_change_possible(rt,st)
    @test !is_change_possible(st,t)
    @test !is_change_possible(st,rt)

    st2 = SkeletonTriangulation(bt)
    test_triangulation(st2)
    st2_normal = get_normal_vector(st2)
    @test is_change_possible(rt,st2)
    @test is_change_possible(bt,st2)

    cell_to_inout = fill(true,num_cells(model1))
    cell_to_inout[1:10] .= false
    it  = InterfaceTriangulation(model1,cell_to_inout)
    @test isa(it,AdaptedTriangulation)
    test_triangulation(it)
    @test is_change_possible(t,it)
    @test is_change_possible(rt,it)
    @test !is_change_possible(it,t)
    @test !is_change_possible(it,rt)
  end
end

end