module PatchTriangulationsTests

using Test
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.Arrays
using Gridap.CellData

using Gridap.Geometry: get_faces
using Gridap.Geometry: PatchTopology, PatchTriangulation
using Gridap.Geometry: PatchBoundaryTriangulation, PatchSkeletonTriangulation
using Gridap.Geometry: InterfacePatchTopology
using Gridap.Geometry: num_patches, get_patch_cells, PolytopalDiscreteModel
using Gridap.Io

model = CartesianDiscreteModel((0,1,0,1),(4,4))
topo = get_grid_topology(model)

# Cell-centered PatchTopology: 1 cell, 4 facets, 4 nodes

ptopo = PatchTopology(model)

ptopo_view = view(ptopo,collect(1:4))
ptopo_rest, _ = Geometry.restrict(ptopo,collect(1:4))
@test all(c -> c ∈ [1,2,3,4], Geometry.get_patch_cells(ptopo_view).data)
@test map(length, Geometry.get_patch_cells(ptopo_view)) == map(length, Geometry.get_patch_cells(ptopo_rest))

patch_cells = Geometry.get_patch_cells(ptopo)
@test all(x -> length(x) == 1, patch_cells)

patch_facets = Geometry.get_patch_facets(ptopo)
@test all(x -> length(x) == 4, patch_facets)

patch_nodes = Geometry.get_patch_faces(ptopo,0)
@test all(x -> length(x) == 4, patch_nodes)

@test Geometry.is_partition(ptopo)
@test Geometry.is_disjoint(ptopo)
@test Geometry.is_cover(ptopo)

ptopo_view = view(ptopo,IdentityVector(num_cells(model)))
ptopo_rest, _ = Geometry.restrict(ptopo,IdentityVector(num_cells(model)))
@test Geometry.get_patch_cells(ptopo_view) == Geometry.get_patch_cells(ptopo_rest) == Geometry.get_patch_cells(ptopo)

Ωp = PatchTriangulation(model,ptopo)
Geometry.test_triangulation(Ωp)

Γp = PatchBoundaryTriangulation(model,ptopo)
Geometry.test_triangulation(Γp)

Γp2 = Geometry.PatchBoundaryTriangulation(model,ptopo;tags="boundary")
Geometry.test_triangulation(Γp2)
@test Geometry.is_change_possible(Γp,Γp2)
@test Geometry.is_change_possible(Γp2,Γp)

fh = CellField(1,Γp)
fh2 = CellData.change_domain(fh,Γp2,ReferenceDomain())
fh3 = CellData.change_domain(fh2,Γp,ReferenceDomain())
@test abs(sum(∫(fh-fh3)*Measure(Γp2,2))) < 1e-12

gh =  CellField(x -> x[1],Γp)
gh2 = CellData.change_domain(gh,Γp2,PhysicalDomain())
gh3 = CellData.change_domain(gh2,Γp,PhysicalDomain())
@test abs(sum(∫(gh-gh3)*Measure(Γp2,2))) < 1e-12

# Vertex-centric PatchTopology: All cells around each vertex

ptopo = PatchTopology(ReferenceFE{0},model)
@test isa(ptopo.metadata,Geometry.StarPatchMetadata)

ptopo_view = view(ptopo,collect(1:4))
ptopo_rest, _ = Geometry.restrict(ptopo,collect(1:4))
@test all(c -> c ∈ [1,2,3,4], Geometry.get_patch_cells(ptopo_view).data)
@test map(length, Geometry.get_patch_cells(ptopo_view)) == map(length, Geometry.get_patch_cells(ptopo_rest))
ptopo_rest, _ = Geometry.restrict(ptopo,collect(1:4);remove_empty_patches=true)

Ωp = PatchTriangulation(model,ptopo)
Γp = PatchBoundaryTriangulation(model,ptopo)
Λp = PatchSkeletonTriangulation(model,ptopo)

ptopo_view = view(ptopo,collect(1:4))
Ω = Triangulation(model,collect(1:4))
Γp1 = PatchBoundaryTriangulation(Ω,ptopo)
Γp2 = PatchBoundaryTriangulation(model,ptopo_view)
@test Γp1.trian.trian.tface_to_mface == Γp2.trian.trian.tface_to_mface
@test Geometry.get_patch_faces(Γp1) == Geometry.get_patch_faces(Γp2)
Λp1 = PatchSkeletonTriangulation(Ω,ptopo)
Λp2 = PatchSkeletonTriangulation(model,ptopo_view)
@test Λp1.trian.plus.trian.tface_to_mface == Λp2.trian.plus.trian.tface_to_mface
@test Geometry.get_patch_faces(Λp1) == Geometry.get_patch_faces(Λp2)

# Interface patch topology

model = CartesianDiscreteModel((0,1,0,1),(2,2))
topo = get_grid_topology(model)

itopo = InterfacePatchTopology(model)
@test Geometry.get_patch_cells(itopo) == [[1],[1],[1],[2],[2],[2],[3],[3],[3],[4],[4],[4]]
@test Geometry.get_patch_faces(itopo,1) == [[1,3],[2],[4],[4],[5,7],[6],[2],[8,9],[10],[6],[10],[11,12]]
Γp = Geometry.PatchBoundaryTriangulation(model,itopo)

itopo = InterfacePatchTopology(
  PatchTopology(topo, Table([[1,2],[3,4]]))
)
@test Geometry.get_patch_cells(itopo) == [[1,2],[1,2],[3,4],[3,4]]
@test Geometry.get_patch_faces(itopo,1) == [[1,3,5,7],[2,6],[2,6],[8,9,11,12]]
Γp = Geometry.PatchBoundaryTriangulation(model,itopo)

# Bounding box grids

model = CartesianDiscreteModel((0,1,0,1),(2,2))
topo = get_grid_topology(model)
patch_cells = Table([[1,2],[3,4]])

ptopo = PatchTopology(topo, patch_cells)
cell_grid = Geometry.bounding_box_grid(model, ptopo, 2)
@test num_cells(cell_grid) == num_patches(ptopo)
@test get_cell_coordinates(cell_grid)[1] == Point.([(0.0, 0.0), (1.0, 0.0), (0.0, 0.5), (1.0, 0.5)])
@test get_cell_coordinates(cell_grid)[2] == Point.([(0.0, 0.5), (1.0, 0.5), (0.0, 1.0), (1.0, 1.0)])

itopo = InterfacePatchTopology(PatchTopology(topo, patch_cells))
Γp = Geometry.PatchBoundaryTriangulation(model,itopo)
face_grid = Geometry.bounding_box_grid(model, Γp; δmin=0.1)
@test num_cells(face_grid) == num_patches(itopo)

# IO

model_io = CartesianDiscreteModel((0,1,0,1),(3,3))
ptopo = PatchTopology(model_io)
dict = to_dict(ptopo)
@test dict[:topo_type] == "UnstructuredGridTopology"
@test haskey(dict,:patch_cells)
ptopo_r = from_dict(PatchTopology,dict)
@test num_patches(ptopo_r) == num_patches(ptopo)
@test num_cells(ptopo_r) == num_cells(ptopo)
pc_orig = get_patch_cells(ptopo)
pc_rest = get_patch_cells(ptopo_r)
@test pc_orig.data == pc_rest.data
@test pc_orig.ptrs == pc_rest.ptrs

pmodel_io = PolytopalDiscreteModel(model_io)
ptopo_p = PatchTopology(pmodel_io)
dict_p = to_dict(ptopo_p)
@test dict_p[:topo_type] == "PolytopalGridTopology"
@test haskey(dict_p,:patch_cells)
ptopo_pr = from_dict(PatchTopology,dict_p)
@test num_patches(ptopo_pr) == num_patches(ptopo_p)
@test num_cells(ptopo_pr) == num_cells(ptopo_p)
pc_orig_p = get_patch_cells(ptopo_p)
pc_rest_p = get_patch_cells(ptopo_pr)
@test pc_orig_p.data == pc_rest_p.data
@test pc_orig_p.ptrs == pc_rest_p.ptrs

# Patches of lower-dimensional faces: one patch per d-face, with its (d-1)-faces as boundary.
# The normal of the boundary is the outward conormal, tangent to the d-face.

function test_face_patches(model,d)
  topo = get_grid_topology(model)
  ptopo = PatchTopology(ReferenceFE{d},ReferenceFE{d},model)
  @test num_cell_dims(ptopo) == d
  @test num_patches(ptopo) == num_faces(topo,d)
  # The sub-topology keeps the numbering of the model
  for n in 0:d-1
    @test get_faces(ptopo.topo,d,n) == get_faces(topo,d,n)
    @test Geometry.get_patch_faces(ptopo,n) == map(sort,get_faces(topo,d,n))
  end

  Fp = PatchTriangulation(model,ptopo)
  Ep = PatchBoundaryTriangulation(model,ptopo)
  @test num_cell_dims(Fp) == d
  @test num_cell_dims(Ep) == d-1
  @test num_cells(Ep) == length(get_faces(topo,d,d-1).data)
  @test isnothing(get_glue(Ep,Val(num_cell_dims(model)))) == (d != num_cell_dims(model))

  # Divergence theorem on each patch: ∫_∂F x⋅n = d |F|
  nE = get_normal_vector(Ep)
  xE = CellField(identity,Ep)
  flux = get_array(∫(xE⋅nE)*Measure(Ep,2))
  meas = get_array(∫(1)*Measure(Fp,2))
  patch_to_tfaces = Ep.glue.patch_to_tfaces
  @test all(p -> sum(flux[patch_to_tfaces[p]]) ≈ d*meas[p], 1:num_patches(ptopo))
end

model2 = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(3,3)))
model3 = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2)))
for model in (model2, Geometry.voronoi(simplexify(model2)))
  test_face_patches(model,1)
  test_face_patches(model,2)
end
for model in (model3, simplexify(model3), PolytopalDiscreteModel(model3))
  test_face_patches(model,1)
  test_face_patches(model,2)
  test_face_patches(model,3)
end
# NOTE: d = 3 is not tested here. The polytopes of a polytopal model created from simplices
# have half of their facet normals pointing inwards (also seen in `BoundaryTriangulation`).
model = PolytopalDiscreteModel(simplexify(model3))
test_face_patches(model,1)
test_face_patches(model,2)

@test PatchTopology(ReferenceFE{0},model3).topo === get_grid_topology(model3)
@test get_patch_cells(PatchTopology(ReferenceFE{2},ReferenceFE{0},model3)) == get_faces(get_grid_topology(model3),0,2)

# Faces of the patch cells of any codimension, glued to the cells

function test_patch_faces(model,df)
  D = num_cell_dims(model)
  topo = get_grid_topology(model)
  ptopo = PatchTopology(ReferenceFE{D},ReferenceFE{D},model)
  Λ = PatchFaceTriangulation(ReferenceFE{df},model,ptopo)
  @test num_cell_dims(Λ) == df
  @test num_cells(Λ) == length(get_faces(topo,D,df).data)
  @test get_glue(Λ,Val(D)).tface_to_mface == Λ.glue.tface_to_patch
  @test get_glue(Λ,Val(df)).tface_to_mface == Geometry.get_patch_faces(ptopo,df).data

  # Quantities on the cells, restricted to the faces
  Ω = Triangulation(model)
  dΛ = Measure(Λ,2)
  meas = get_array(∫(1)dΛ)
  id = CellField(collect(Float64,1:num_cells(Ω)),Ω)
  @test get_array(∫(id)dΛ) ≈ Λ.glue.tface_to_patch .* meas
  x = CellData.GenericCellField(get_cell_map(Ω),Ω,ReferenceDomain())
  @test get_array(∫(x⋅x)dΛ) ≈ get_array(∫(CellField(y -> y⋅y,Λ))dΛ)
end

for model in (model2, Geometry.voronoi(simplexify(model2)))
  test_patch_faces(model,0)
  test_patch_faces(model,1)
end
for model in (model3, simplexify(model3), PolytopalDiscreteModel(model3))
  test_patch_faces(model,0)
  test_patch_faces(model,1)
  test_patch_faces(model,2)
end

end # module
