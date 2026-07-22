module LoopPatchVerticesMapTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Visualization
using GridapSubdivisionSurfaces.FESpaces: LoopPatchVerticesMap

const tmp_dir = joinpath(@__DIR__, "..", "tmp")
mkpath(tmp_dir)
tmp_path(name) = joinpath(tmp_dir, name)

nverts_core(N) = 1 + N + 2N
ncells_core(N) = 4N

"""
    isolated_ev_mesh(N)

Builds a small triangulated mesh with a single isolated extraordinary vertex
(EV) valence `N`, within a regular 2-ring neighborhood.

Returns `(topo, central_triangle)` with `central_triangle = [1,2,3]`.
"""
function isolated_ev_mesh(N)
  ring1 = collect(2:N+1)          # N vertices
  ring2 = collect(N+2:N+1+2N)     # 2N vertices

  conn = Vector{Int32}[]

  # EV 1-neighborhood / first ring: EV(1) -- ring1
  for i in 1:N
    a = ring1[i]; b = ring1[mod1(i+1, N)]
    push!(conn, Int32[1, a, b])
  end

  # EV second ring: ring1 -- ring2 (gives ring1 vertices valence 6)
  for i in 1:N
    a = ring1[i]; b = ring1[mod1(i+1, N)]
    r2a = ring2[mod1(2i-1, 2N)]; r2b = ring2[mod1(2i, 2N)]; r2c = ring2[mod1(2i+1, 2N)]
    push!(conn, Int32[a, r2a, r2b])
    push!(conn, Int32[a, r2b, b])
    push!(conn, Int32[b, r2b, r2c])
  end

  # close every ring2 vertex to valence 6 with dummy neighboors
  next_dummy = nverts_core(N) + 1
  for k in 1:2N
    v = ring2[k]
    ndummies = isodd(k) ? 3 : 2
    dummies = collect(next_dummy:next_dummy+ndummies-1)
    next_dummy += ndummies
    for j in 1:ndummies-1
      push!(conn, Int32[v, dummies[j], dummies[j+1]])
    end
  end
  nverts = next_dummy - 1

  cell_vertices = Table(conn)
  reffes = [LagrangianRefFE(Float64, TRI, 1)]
  cell_types = ones(Int8, length(conn))

  # real 2D layout for the core (EV at the origin, ring1 at radius 1, ring2
  # at radius 2), for visualization (see export_patch_vtk)
  coords = fill(Point(0.0, 0.0), nverts)
  for i in 1:N
    θ = 2π*(i-1)/N
    coords[ring1[i]] = Point(cos(θ), sin(θ))
  end
  for k in 1:2N
    θ = π*(k-1)/N
    coords[ring2[k]] = Point(2*cos(θ), 2*sin(θ))
  end

  grid = UnstructuredGrid(coords, cell_vertices, reffes, cell_types)
  topo = UnstructuredGridTopology(grid)

  topo, Int32[1, 2, 3]
end

# exports the EV patch as a triangulation
function export_patch_vtk(N, filebase)
  topo, _ = isolated_ev_mesh(N)
  coords = get_vertex_coordinates(topo)[1:nverts_core(N)]
  core_cell_vertices = get_cell_vertices(topo)[1:ncells_core(N)]

  reffes = [LagrangianRefFE(Float64, TRI, 1)]
  cell_types = ones(Int8, ncells_core(N))
  grid = UnstructuredGrid(coords, core_cell_vertices, reffes, cell_types)
  core_topo = UnstructuredGridTopology(grid)
  model = UnstructuredDiscreteModel(grid, core_topo, FaceLabeling(core_topo))
  trian = Triangulation(model)

  writevtk(trian, filebase)
end

graph_neighbors(g, v) = @view g.rowval[g.colptr[v]:g.colptr[v+1]-1]
valence(g, v) = g.colptr[v+1] - g.colptr[v]
adjacent(g, u, v) = v ∈ graph_neighbors(g, u)

@testset "LoopPatchVerticesMap, isolated EV of valence $N" for N in (3, 4, 5, 6, 7)
  topo, central_triangle = isolated_ev_mesh(N)
  g = compute_graph(topo, 0, 1; self_loops=false)

  @test valence(g, 1) == N
  @test all(v -> valence(g, v) == 6, 2:N+1)

  loop_patch_map = LoopPatchVerticesMap(topo)
  cache = return_cache(loop_patch_map, central_triangle)
  patch_vertices = collect(evaluate!(cache, loop_patch_map, central_triangle))

  @test allunique(patch_vertices)

  if N == 6
    # regular patch: unaffected, existing code path
    @test length(patch_vertices) == 12
    v4, v7, v8 = patch_vertices[4], patch_vertices[7], patch_vertices[8]
    @test v4 == 1
    @test Set((v7, v8)) == Set((2, 3))
    @test adjacent(g, v4, v7)
    @test adjacent(g, v4, v8)
  else
    # irregular patch: K = N+6 vertices, Stam (1998) Fig. 2 numbering
    @test length(patch_vertices) == N+6

    ev, v2, vN1 = patch_vertices[1], patch_vertices[2], patch_vertices[N+1]
    @test ev == 1
    @test Set((v2, vN1)) == Set((2, 3))

    # ring 2,…,N+1 is a cyclic path in the mesh graph
    for i in 2:N
      @test adjacent(g, patch_vertices[i], patch_vertices[i+1])
    end
    @test adjacent(g, v2, vN1)
    @test all(v -> adjacent(g, ev, v), patch_vertices[2:N+1])

    v3, vN = patch_vertices[3], patch_vertices[N]
    vN2, vN3, vN4, vN5, vN6 = patch_vertices[N+2:N+6]
    @test adjacent(g, vN2, v2) && adjacent(g, vN2, vN1)
    @test adjacent(g, vN4, v2) && adjacent(g, vN4, v3)
    @test adjacent(g, vN6, vN1) && adjacent(g, vN6, vN)
    @test adjacent(g, vN3, v2) && adjacent(g, vN3, vN2)
    @test adjacent(g, vN5, vN1) && adjacent(g, vN5, vN2)
  end
end

end # module
