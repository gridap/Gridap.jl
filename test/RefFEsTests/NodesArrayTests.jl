module NodesArrayTests

using Gridap, Test

p = Polytope(HEX_AXIS)
test_hex_polytopes(p)

p = Polytope(HEX_AXIS, HEX_AXIS)
test_hex_polytopes(p)

p = Polytope(HEX_AXIS, HEX_AXIS, HEX_AXIS)
test_hex_polytopes(p)

p = Polytope(HEX_AXIS,TET_AXIS)
# _high_order_test(p,Float64,1)
T = Float64
horeffe = Gridap.RefFEs._high_order_lagrangian_reffe(p, T, 2)
@test length(horeffe.shfbasis) == 6
horeffe = Gridap.RefFEs._high_order_lagrangian_reffe(p, T, 3)
@test length(horeffe.shfbasis) == 10

p = Polytope(HEX_AXIS,TET_AXIS,TET_AXIS)
# _high_order_test(p,Float64,1)
T = Float64
horeffe = Gridap.RefFEs._high_order_lagrangian_reffe(p, T, 2)
@test length(horeffe.shfbasis) == 10
horeffe = Gridap.RefFEs._high_order_lagrangian_reffe(p, T, 3)
@test length(horeffe.shfbasis) == 20

# Problems with grid for new machinery
p = Polytope(HEX_AXIS,TET_AXIS,TET_AXIS)
LagrangianRefFE(Float64,p,3)
grid = Grid(p,2)
# writevtk(grid,"grid")

function _high_order_test_scalar(p,T,order)
  horeffe = Gridap.RefFEs._high_order_lagrangian_reffe(p, T, order)
  if (Gridap.RefFEs._is_hex(p)) @test length(horeffe.shfbasis) == (order+1)^(length(p.extrusion)) end
  # if (Gridap.RefFEs._is_tet(p)) @test length(horeffe.shfbasis) == sum(1:order+1) end
end

function _high_order_test_vector(p,T,order)
  horeffe = Gridap.RefFEs._high_order_lagrangian_reffe(p, T, order)
  if (Gridap.RefFEs._is_hex(p)) @test length(horeffe.shfbasis) == length(p.extrusion)*(order+1)^(length(p.extrusion)) end
end

function test_hex_polytopes(p)
  _high_order_test_scalar(p,Float64,1)
  _high_order_test_scalar(p,Float64,2)
  _high_order_test_scalar(p,Float64,3)
  _high_order_test_scalar(p,Float64,4)
  T = VectorValue{length(p.extrusion),Float64}
  _high_order_test_vector(p,T,1)
  _high_order_test_vector(p,T,2)
  _high_order_test_vector(p,T,3)
  _high_order_test_vector(p,T,4)
end

end #module
