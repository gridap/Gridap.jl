module PatchTriangulationsTests

using Test
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.Arrays
using Gridap.CellData

using Gridap.Geometry: get_faces
using Gridap.Geometry: PatchTopology, PatchTriangulation
using Gridap.Geometry: PatchBoundaryTriangulation, PatchSkeletonTriangulation
using Gridap.Geometry: InterfacePatchTopology

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

end