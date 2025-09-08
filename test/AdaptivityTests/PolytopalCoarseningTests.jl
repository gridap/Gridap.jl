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

cmodel = Adaptivity.coarsen(model,ptopo)
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

cmodel = Adaptivity.coarsen(model,ptopo)
writevtk(Triangulation(ReferenceFE{1},cmodel),joinpath(outdir,"cmodel");append=false)

for (i,p) in enumerate(get_polytopes(cmodel))
  writevtk(p,joinpath(outdir,"p_$i");append=false)
end

end