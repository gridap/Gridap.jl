module MacroGridTopologiesTests

using Test
using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Geometry: get_faces, get_fine_to_coarse_faces, get_coarse_to_fine_faces
using Gridap.Geometry: get_coarse_vertex_to_vertex
using Gridap.Geometry: get_face_coordinates, get_vertex_coordinates, get_isboundary_face
using Gridap.Geometry: get_patch_cells, MacroStarPatchMetadata
using Gridap.FESpaces

# Blocks of b^D cells of a Cartesian mesh with n^D cells
function block_partition(n,b,D)
  nb = cld(n,b)
  cell_to_ccell = Int32[]
  for I in CartesianIndices(ntuple(_->n,D))
    J = map(i->(i-1)÷b,Tuple(I))
    push!(cell_to_ccell,1+sum(J[k]*nb^(k-1) for k in 1:D))
  end
  return cell_to_ccell
end

# One tag per side of the unit box on the boundary facets (0 elsewhere)
function side_tags(topo)
  D = num_cell_dims(topo)
  isboundary = get_isboundary_face(topo,D-1)
  map(get_face_coordinates(topo,D-1),isboundary) do xs, isb
    isb || return Int32(0)
    x = sum(xs)/length(xs)
    k = findfirst(k -> x[k] < 1.0e-10 || x[k] > 1-1.0e-10, 1:D)
    Int32(2k - (x[k] < 0.5))
  end
end

sorted_rows(t) = [sort(collect(r)) for r in t]

for D in (2,3)
  n = 4
  model = CartesianDiscreteModel(Tuple(repeat([0,1],D)),ntuple(_->n,D))
  topo = get_grid_topology(model)
  tags = side_tags(topo)

  # Every cell is a coarse cell, every boundary facet its own tag: same topology
  nfacets = num_faces(topo,D-1)
  mt = MacroGridTopology(topo;facet_to_tag=collect(Int32,1:nfacets))
  @test get_grid_topology(mt) === topo
  @test all(isassigned(mt.n_m_to_nface_to_mfaces,i) == false for i in eachindex(mt.n_m_to_nface_to_mfaces))
  for d in 0:D
    @test num_faces(mt,d) == num_faces(topo,d)
    cdims, cfaces = get_fine_to_coarse_faces(mt,d)
    @test all(==(d),cdims) && cfaces == 1:num_faces(topo,d)
  end
  for n in 0:D, m in 0:D
    @test sorted_rows(get_faces(mt,n,m)) == sorted_rows(get_faces(topo,n,m))
  end

  # Blocks of 2^D cells, one tag per side: a Cartesian mesh with 2^D cells
  cell_to_ccell = block_partition(n,2,D)
  mt = MacroGridTopology(topo,cell_to_ccell;facet_to_tag=tags)
  ctopo = get_grid_topology(CartesianDiscreteModel(Tuple(repeat([0,1],D)),ntuple(_->2,D)))
  @test [num_faces(mt,d) for d in 0:D] == [num_faces(ctopo,d) for d in 0:D]
  @test sort(get_vertex_coordinates(mt)) == sort(get_vertex_coordinates(ctopo))
  @test count(get_isboundary_face(mt,D-1)) == count(get_isboundary_face(ctopo,D-1))
  @test length(get_coarse_to_fine_faces(mt,D,D)) == 2^D
  @test all(length(cells) == 2^D for cells in get_coarse_to_fine_faces(mt,D,D))
  @test all(length(facets) == 2^(D-1) for facets in get_coarse_to_fine_faces(mt,D-1,D-1))
  @test all(length(cvertices) == 2^D for cvertices in get_faces(mt,D,0))
  @test length(get_coarse_vertex_to_vertex(mt)) == num_faces(mt,0)
  @test sum((-1)^d*num_faces(mt,d) for d in 0:D) == 1

  # Default tags: the boundary of each coarse cell is a single coarse facet (corners are lost)
  mt = MacroGridTopology(topo,cell_to_ccell)
  if D == 2
    @test [num_faces(mt,d) for d in 0:D] == [5,8,4]
  end
  @test sum((-1)^d*num_faces(mt,d) for d in 0:D) == 1
end

# Patch topologies
for D in (2,3)
  n = 4
  model = CartesianDiscreteModel(Tuple(repeat([0,1],D)),ntuple(_->n,D))
  topo = get_grid_topology(model)
  cell_to_ccell = block_partition(n,2,D)
  mt = MacroGridTopology(topo,cell_to_ccell;facet_to_tag=side_tags(topo))

  # Coarse cells
  ptopo = PatchTopology(mt)
  @test ptopo.metadata isa MacroStarPatchMetadata
  @test sorted_rows(get_patch_cells(ptopo)) == sorted_rows(get_coarse_to_fine_faces(mt,D,D))

  # One patch per coarse facet, with its fine facets
  ptopo = PatchTopology(ReferenceFE{D-1},ReferenceFE{D-1},mt)
  @test num_cells(ptopo.topo) == num_faces(topo,D-1)
  @test sorted_rows(get_patch_cells(ptopo)) == sorted_rows(get_coarse_to_fine_faces(mt,D-1,D-1))

  # Coarse vertex stars: fine cells of the coarse cells around each coarse vertex
  ptopo = PatchTopology(ReferenceFE{0},mt)
  @test Geometry.num_patches(ptopo) == num_faces(mt,0)
  @test all(length(cells) == 2^D*length(cc) for (cells,cc) in zip(get_patch_cells(ptopo),get_faces(mt,0,D)))
  ptopo = PatchTopology(ReferenceFE{0},mt;patch_roots=[1,3])
  @test ptopo.metadata.patch_roots == [1,3]

  # Every cell a coarse cell, every boundary facet its own tag: :star on the coarse vertex stars
  # is the :star of the vertex stars of the model
  mt = MacroGridTopology(topo;facet_to_tag=collect(Int32,1:num_faces(topo,D-1)))
  V = FESpace(model,ReferenceFE(lagrangian,Float64,2);dirichlet_tags="boundary")
  rows_macro = FESpaces.get_patch_assembly_ids(Val(:star),V,PatchTopology(ReferenceFE{0},mt))
  rows = FESpaces.get_patch_assembly_ids(Val(:star),V,PatchTopology(ReferenceFE{0},model))
  @test sorted_rows(rows_macro) == sorted_rows(rows)
end

end # module
