module NonOrientedMomentBasedTests

# Raviart–Thomas, Nédélec and BDM on meshes whose cells list the vertices of
# shared faces in different orders: the DoF numbering permutes the face DoFs
# and the change of basis applies their signs, so that a field of the space
# is interpolated exactly.

using Test
using Random
using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.FESpaces

# 2×2 QUAD mesh on [0,2]² with one cell flipped: the edge between cells 1 and
# 2 has pindex 2 from cell 2
function unoriented_quad_model()
  x = Point{2,Float64}[(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)]
  cells = Table(Int32[1,2,4,5, 5,6,2,3, 4,5,7,8, 5,6,8,9], Int32[1,5,9,13,17])
  reffe = ReferenceFE(QUAD, lagrangian, Float64, 1)
  UnstructuredDiscreteModel(UnstructuredGrid(x, cells, [reffe], fill(Int8(1), 4), NonOriented()))
end

# 2×2×2 HEX mesh on [0,2]³ with one cell's y and z axes swapped: two of its
# faces have pindex 2 from the neighbours, and no edge is reversed
function unoriented_hex_model()
  nidx(i,j,k) = Int32(i + 3j + 9k + 1)
  x = Point{3,Float64}[(i,j,k) for k in 0:2 for j in 0:2 for i in 0:2]
  hex(i,j,k) = Int32[nidx(i,j,k), nidx(i+1,j,k), nidx(i,j+1,k), nidx(i+1,j+1,k),
                     nidx(i,j,k+1), nidx(i+1,j,k+1), nidx(i,j+1,k+1), nidx(i+1,j+1,k+1)]
  swapped = Int32[4,13,5,14, 7,16,8,17]
  cells = vcat(hex(0,0,0), hex(1,0,0), swapped, hex(1,1,0), hex(0,0,1), hex(1,0,1), hex(0,1,1), hex(1,1,1))
  reffe = ReferenceFE(HEX, lagrangian, Float64, 1)
  UnstructuredDiscreteModel(UnstructuredGrid(x, Table(cells, Int32[1 + 8i for i in 0:8]), [reffe], fill(Int8(1), 8), NonOriented()))
end

# the same mesh with the vertices of each cell relabelled by a symmetry of its
# polytope, in a scrambled order
function relabel_cells(model)
  grid = get_grid(model)
  syms = get_vertex_permutations(only(get_polytopes(model)))
  perms = syms[randperm(MersenneTwister(0), length(syms))]
  cell_nodes = get_cell_node_ids(grid)
  relabelled = [collect(cell_nodes[i])[perms[mod1(i, length(perms))]] for i in 1:length(cell_nodes)]
  x = vec(collect(get_node_coordinates(grid)))
  UnstructuredDiscreteModel(UnstructuredGrid(
    x, Table(relabelled), collect(get_reffes(grid)), collect(get_cell_type(grid)), NonOriented()))
end

function interpolation_error(model, reffe, u)
  V = FESpace(model, reffe)
  uh = interpolate(u, V)
  dΩ = Measure(Triangulation(model), 2*(get_order(reffe) + 2))
  sqrt(sum(∫((uh - u)⋅(uh - u))dΩ))
end

# a field of every space of order ≥ 1, and of order 0
linear(D) = x -> VectorValue(ntuple(i -> 1.0 + 2.0*x[i] - 0.5*sum(x), D))
constant(D) = x -> VectorValue(ntuple(i -> 1.0 - 0.5*i, D))

tol = 1e-10

for (p, model) in ((QUAD, unoriented_quad_model()), (HEX, unoriented_hex_model()))
  D = num_dims(p)
  for order in 0:2, name in (raviart_thomas, nedelec)
    reffe = ReferenceFE(p, name, order; change_dof=false)
    u = order == 0 ? constant(D) : linear(D)
    @test interpolation_error(model, reffe, u) < tol
  end
end

# simplices: every vertex permutation of a cell
for D in (2, 3)
  domain = ntuple(i -> isodd(i) ? 0 : 1, 2*D)
  model = relabel_cells(simplexify(CartesianDiscreteModel(domain, ntuple(_ -> 2, D))))
  p = D == 2 ? TRI : TET
  for order in 0:2
    u = order == 0 ? constant(D) : linear(D)
    @test interpolation_error(model, ReferenceFE(p, raviart_thomas, order; change_dof=false), u) < tol
    order > 0 && @test interpolation_error(model, ReferenceFE(p, bdm, order; change_dof=false), u) < tol
  end
  # Nédélec's triangular face weights have no permutation tables yet
  nd_orders = D == 2 ? (0:2) : (0:0)
  for order in nd_orders
    u = order == 0 ? constant(D) : linear(D)
    @test interpolation_error(model, ReferenceFE(p, nedelec, order; change_dof=false), u) < tol
  end
end

end # module
