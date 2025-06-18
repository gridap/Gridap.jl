module PatchTriangulationsTests

using Test
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.Arrays
using Gridap.CellData

model = CartesianDiscreteModel((0,1,0,1),(4,4))

ptopo = Geometry.PatchTopology(model)

patch_cells = Geometry.get_patch_cells(ptopo)
patch_facets = Geometry.get_patch_facets(ptopo)
patch_nodes = Geometry.get_patch_faces(ptopo,0)

Ωp = Geometry.PatchTriangulation(model,ptopo)
Geometry.test_triangulation(Ωp)

Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)
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

end