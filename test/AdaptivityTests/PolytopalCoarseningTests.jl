module PolytopalCoarseningTests

using Test
using Gridap
using Gridap.ReferenceFEs, Gridap.Adaptivity, Gridap.Geometry, Gridap.Arrays

function select_facet(p,lf,cperms)
  D = num_dims(p)
  f = get_faces(p,D-1,0)[lf]
  fperms = get_vertex_permutations(Polytope{D-1}(p,lf))
  pid = cperms[lf + get_offset(p,D-1)]
  return f[fperms[pid]]
end

outdir = "tmp/"
mkpath(outdir)

for D in 2:3
  n = 2
  nc = (n,fill(1,D-1)...)
  domain = (0,n,repeat([0,1],D-1)...)
  model = Geometry.PolytopalDiscreteModel(CartesianDiscreteModel(domain,nc))

  topo = get_grid_topology(model)
  ptopo = Geometry.PatchTopology(topo,Table([collect(1:n)]))

  cell_polys = get_polytopes(topo)
  cell_perms = get_cell_permutations(topo)
  cell_to_facets = get_faces(topo,D,D-1)
  facet_to_cells = get_faces(topo,D-1,D)

  facet = findfirst(x -> length(x) == 2, facet_to_cells)
  c1, c2 = facet_to_cells[facet]
  p1, p2 = cell_polys[c1], cell_polys[c2]
  lf1, lf2 = findfirst(isequal(facet), cell_to_facets[c1]), findfirst(isequal(facet), cell_to_facets[c2])
  f1, f2 = select_facet(p1,lf1,cell_perms[c1]), select_facet(p2,lf2,cell_perms[c2])

  p = ReferenceFEs.merge_polytopes(p1, p2, f1, f2)
  q, _ = ReferenceFEs.polytope_from_faces(D,get_vertex_coordinates(p),get_faces(p,D-1,0))

  writevtk(p1,joinpath(outdir,"p1");append=false)
  writevtk(p2,joinpath(outdir,"p2");append=false)
  writevtk(p,joinpath(outdir,"p");append=false)
  writevtk(q,joinpath(outdir,"q");append=false)
end

############################################################################################

D = 2
n = 4
nc = (n,fill(1,D-1)...)
domain = (0,n,repeat([0,1],D-1)...)
model = Geometry.voronoi(simplexify(CartesianDiscreteModel(domain,nc)))
writevtk(Triangulation(model),joinpath(outdir,"model");append=false)

patch_cells = Table([
  [1],[2,3,6,7],[4,8],[5,9,10]
])
topo = get_grid_topology(model)
ptopo = Geometry.PatchTopology(topo,patch_cells)

cmodel, glue = Adaptivity.coarsen(model,ptopo;return_glue=true)
writevtk(Triangulation(cmodel),joinpath(outdir,"cmodel");append=false)

############################################################################################

D = 3
nc = (3,3,3)
domain = (0,1,0,1,0,1)
model = Geometry.PolytopalDiscreteModel(CartesianDiscreteModel(domain,nc))
writevtk(Triangulation(ReferenceFE{1},model),joinpath(outdir,"model");append=false)

patch_cells = Table([
  [1],[2,3,4,5,6,7,8,9],[10,11,13,14,19,20,22,23],[12,15,16,17,18],[21,24],[25,26],[27]
])
topo = get_grid_topology(model)
ptopo = Geometry.PatchTopology(topo,patch_cells)

cmodel, glue = Adaptivity.coarsen(model,ptopo;return_glue=true)
writevtk(Triangulation(ReferenceFE{1},cmodel),joinpath(outdir,"cmodel");append=false)

for (i,p) in enumerate(get_polytopes(cmodel))
  writevtk(p,joinpath(outdir,"p_$i");append=false)
end

############################################################################################
# Glue composition

model1 = Geometry.voronoi(Geometry.simplexify(CartesianDiscreteModel((0,1,0,1), (3,3))))

pcells1 = Table([
  LinearIndices((4,4))[1:2,1:2],
  LinearIndices((4,4))[3:4,1:2],
  LinearIndices((4,4))[1:2,3:4],
  LinearIndices((4,4))[3:4,3:4],
])
ptopo1 = Geometry.PatchTopology(get_grid_topology(model1), pcells1)
model2, glue12 = Adaptivity.coarsen(model1, ptopo1; return_glue=true)

pcells2 = Table([[1,2],[3,4]])
ptopo2 = Geometry.PatchTopology(get_grid_topology(model2), pcells2)
model3, glue23 = Adaptivity.coarsen(model2, ptopo2; return_glue=true)

glue13 = Adaptivity.compose_glues(glue12, glue23)

Dc = 2
topo1 = get_grid_topology(model1)
topo3 = get_grid_topology(model3)
for d in 0:Dc
  n2o_dfaces = glue13.n2o_faces_map[d+1]
  @test length(n2o_dfaces) == num_faces(topo1, d)
  if d == Dc
    @test count(iszero, n2o_dfaces) == 0
  else
    @test count(!iszero, n2o_dfaces) == num_faces(topo3, d)
  end
end

end