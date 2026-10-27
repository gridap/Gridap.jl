module MomentBasedQuadraturesTests

using Test
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.CellData
using Gridap.FESpaces, Gridap.Adaptivity

# A polygon is convex iff all its consecutive edges turn in the same direction
function is_convex(p::Polygon)
  x = get_vertex_coordinates(p)
  n = length(x)
  turns = map(1:n) do i
    a, b, c = x[i], x[mod1(i+1,n)], x[mod1(i+2,n)]
    (b-a)[1]*(c-b)[2] - (b-a)[2]*(c-b)[1]
  end
  all(t -> t >= -1e-12, turns) || all(t -> t <= 1e-12, turns)
end

monomials(D, degree, filter) = [
  x -> prod(ntuple(i -> x[i]^e[i], D))
  for e in Iterators.product(ntuple(i -> 0:degree, D)...) if filter(e, degree)
]

q_filter = (e, o) -> (maximum(e) <= o)
p_filter = (e, o) -> (sum(e) <= o)

function patch_integrals(f, dΩ, patch_cells)
  v = collect(get_array(∫(f)dΩ))
  [sum(v[cells]) for cells in patch_cells]
end

function test_patch_exactness(dΩc, dΩf, patch_cells, fs)
  for f in fs
    vc = collect(get_array(∫(f)dΩc))
    vf = patch_integrals(f, dΩf, patch_cells)
    @test maximum(abs.(vc .- vf)) < 1e-12
  end
end

############################################################################################
# 2D: L-shaped and U-shaped patches
#
#  13 14 15 16       5 6 5 7
#   9 10 11 12  -->  5 5 5 4
#   5  6  7  8       1 3 4 4
#   1  2  3  4       1 1 2 2

model = Geometry.PolytopalDiscreteModel(CartesianDiscreteModel((0,1,0,1),(4,4)))
patch_cells = Table([[1,2,5],[3,4],[6],[7,8,12],[9,10,11,13,15],[14],[16]])
ptopo = Geometry.PatchTopology(get_grid_topology(model),patch_cells)
cmodel = Adaptivity.coarsen(model,ptopo)

polys = get_polytopes(cmodel)
@test !is_convex(polys[1]) # L
@test is_convex(polys[2])
@test !is_convex(polys[4]) # L
@test !is_convex(polys[5]) # U

Ωf = Triangulation(model)
Ωc = Triangulation(cmodel)

for degree in 1:3
  # Fine quadrature exact for Q_degree
  dΩf = Measure(Ωf, 2*degree)

  # Patch measure: gathers the fine quadratures, exact for any integrand
  dΩp = PatchMeasure(Ωc, dΩf, ptopo)
  test_patch_exactness(dΩp, dΩf, patch_cells, [x -> sin(x[1])*exp(x[2])])

  # Moment-based measure, Q space
  dΩm = MomentBasedMeasure(Ωc, dΩf, ptopo, degree)
  quad = CellData.get_cell_quadrature(dΩm)
  @test all(isequal((degree+1)^2), map(length, quad.cell_point))
  test_patch_exactness(dΩm, dΩf, patch_cells, monomials(2, degree, q_filter))
  @test sum(∫(1)dΩm) ≈ 1.0

  # Moment-based measure, P space
  dΩm = MomentBasedMeasure(Ωc, dΩf, ptopo, degree; space=:P)
  quad = CellData.get_cell_quadrature(dΩm)
  @test all(isequal(binomial(degree+2,2)), map(length, quad.cell_point))
  test_patch_exactness(dΩm, dΩf, patch_cells, monomials(2, degree, p_filter))

  # Moment-based measure from the patch measure
  dΩm2 = MomentBasedMeasure(Ωc, dΩp, degree)
  test_patch_exactness(dΩm2, dΩf, patch_cells, monomials(2, degree, q_filter))
end

# Fine measure in the reference domain (non-polytopal fine triangulation)
degree = 2
cart = CartesianDiscreteModel((0,1,0,1),(4,4))
dΩf = Measure(Triangulation(cart), 2*degree)
dΩm = MomentBasedMeasure(Ωc, dΩf, ptopo, degree)
test_patch_exactness(dΩm, dΩf, patch_cells, monomials(2, degree, q_filter))

# Mass matrix of a polytopal FE space: P1 x P1 ⊂ Q2
dΩf = Measure(Ωf, 2*degree)
dΩp = PatchMeasure(Ωc, dΩf, ptopo)
dΩm = MomentBasedMeasure(Ωc, dΩf, ptopo, degree)
V = FESpaces.PolytopalFESpace(cmodel, Float64, 1; space=:P)
Mp = assemble_matrix((u,v) -> ∫(u*v)dΩp, V, V)
Mm = assemble_matrix((u,v) -> ∫(u*v)dΩm, V, V)
@test maximum(abs.(Mp - Mm)) < 1e-12

# Fine triangulation that does not cover the whole fine model (cell 2 removed)
fcells = setdiff(1:16, [2])
Ωf_sub = Triangulation(model, fcells)
dΩf_sub = Measure(Ωf_sub, 2*degree)
dΩm = MomentBasedMeasure(Ωc, dΩf_sub, ptopo, degree)
fcell_to_ftface = Geometry.get_glue(Ωf_sub, Val(2)).mface_to_tface
patch_ftfaces = [filter(>(0), fcell_to_ftface[cells]) for cells in patch_cells]
test_patch_exactness(dΩm, dΩf_sub, patch_ftfaces, monomials(2, degree, q_filter))
@test sum(∫(1)dΩm) ≈ 15/16

############################################################################################
# 2D, embedded-like: the coarse cells are (reference) QUADs of a background mesh, and the
# fine cells are triangles covering an L-shaped part of some of them.
#
# The coarse-to-fine table is indexed by triangulation ids, and the moment fit is done on
# the physical bounding box of each coarse cell.

nc = 4
bgmodel = CartesianDiscreteModel((0,1,0,1),(nc,nc))
fmodel = simplexify(CartesianDiscreteModel((0,1,0,1),(2*nc,2*nc)))

coarse_cell(x) = (floor(Int,x[1]*nc)+1) + nc*floor(Int,x[2]*nc)
is_cut(cell) = (cell-1) % nc == (cell-1) ÷ nc  # cells on the diagonal
# Remove the fine triangles in the top-right quarter of the cut cells
function is_inside(x)
  c = coarse_cell(x)
  !(is_cut(c) && x[1]*nc - floor(x[1]*nc) > 0.5 && x[2]*nc - floor(x[2]*nc) > 0.5)
end

fcell_centroids = lazy_map(xs -> sum(xs)/length(xs), get_cell_coordinates(Triangulation(fmodel)))
Ωf = Triangulation(fmodel, findall(is_inside, collect(fcell_centroids)))
ftface_centroids = lazy_map(xs -> sum(xs)/length(xs), get_cell_coordinates(Ωf))

# Coarse triangulation: only the cut cells (active triangulation ids != model ids)
Ωc = Triangulation(bgmodel, findall(is_cut, 1:nc*nc))
ctface_to_cell = Geometry.get_glue(Ωc, Val(2)).tface_to_mface
cell_to_ctface = Geometry.get_glue(Ωc, Val(2)).mface_to_tface
ctface_to_ftfaces = [Int32[] for _ in 1:num_cells(Ωc)]
for (ftface, x) in enumerate(ftface_centroids)
  ctface = cell_to_ctface[coarse_cell(x)]
  ctface > 0 && push!(ctface_to_ftfaces[ctface], ftface)
end
ctface_to_ftfaces = Table(ctface_to_ftfaces)

for degree in 1:3
  dΩf = Measure(Ωf, 2*degree)

  dΩp = PatchMeasure(Ωc, dΩf, ctface_to_ftfaces)
  test_patch_exactness(dΩp, dΩf, ctface_to_ftfaces, [x -> sin(x[1])*exp(x[2])])
  @test sum(∫(1)dΩp) ≈ num_cells(Ωc) * (3/4) / nc^2

  dΩm = MomentBasedMeasure(Ωc, dΩf, ctface_to_ftfaces, degree)
  test_patch_exactness(dΩm, dΩf, ctface_to_ftfaces, monomials(2, degree, q_filter))

  # Nodes are on the physical cells, not on the reference QUAD
  quad = CellData.get_cell_quadrature(dΩm)
  for (ctface, xs) in enumerate(quad.cell_point)
    pmin, pmax = get_bounding_box(get_cell_coordinates(Ωc)[ctface])
    @test all(x -> all(pmin.data .- 1e-12 .<= x.data .<= pmax.data .+ 1e-12), xs)
  end
end

############################################################################################
# 3D: L-shaped patch

model = Geometry.PolytopalDiscreteModel(CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2)))
patch_cells = Table([[1,2,3],[4],[5,6,7,8]])
ptopo = Geometry.PatchTopology(get_grid_topology(model),patch_cells)
cmodel = Adaptivity.coarsen(model,ptopo)
@test get_measure.(get_polytopes(cmodel)) ≈ [0.375, 0.125, 0.5]

Ωf = Triangulation(model)
Ωc = Triangulation(cmodel)
for degree in 1:2
  dΩf = Measure(Ωf, 3*degree)
  dΩm = MomentBasedMeasure(Ωc, dΩf, ptopo, degree)
  test_patch_exactness(dΩm, dΩf, patch_cells, monomials(3, degree, q_filter))
  dΩm = MomentBasedMeasure(Ωc, dΩf, ptopo, degree; space=:P)
  test_patch_exactness(dΩm, dΩf, patch_cells, monomials(3, degree, p_filter))
end

end # module
