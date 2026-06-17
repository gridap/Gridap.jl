module LoopSurfaceTests

using Test
using Gridap
using GridapSubdivisionSurfaces
using GridapSubdivisionSurfaces.ReferenceFEs

using Gridap.ReferenceFEs

using Gridap.Geometry
using Gridap.Geometry: get_patch_faces
using Gridap.Geometry: get_patch_cells
using Gridap.Geometry: get_patch_to_tfaces
using Gridap.Geometry: get_pface_to_lpface
using Gridap.Geometry: get_pface_to_patch
using Gridap.Geometry: num_patches
using Gridap.Geometry: extend_patches_by_single_layer

"""
    simplex_torus_chart(n=4)

Torus chart grid ([0,1]²) meshed with 2*n² triangles)
"""
function Loop_torus_chart_and_map(n1=10, n2=2*n1)
  model = CartesianDiscreteModel((0,1,0,1), (n1,n2), isperiodic=(true,true))
  model = simplexify(model)
  torus_chart = UnstructuredDiscreteModel(model)

  function torus_map(p)
    x,y = p
    R = 5
    r = 1
    θ = 2π*x
    φ = 2π*y
    Point( (R+r*sin(θ))*cos(φ), (R+r*sin(θ))*sin(φ), r*cos(θ) )
  end

  torus_chart, torus_map
end
torus_chart, torus_map = Loop_torus_chart_and_map(6)
#writevtk(torus_chart , "torus_chart")

loop_torus_model = loop_surface_model(torus_chart, torus_map)
loop_torus_trian = Triangulation(loop_torus_model)
writevtk(loop_torus_trian, "loop_torus"; nsubcells=20)

control_vertices = get_vertex_coordinates(get_grid_topology(loop_torus_model ))
writevtk(control_vertices, "control_vertices")

#T = Float64
#chart_model = torus_chart
#topo = get_grid_topology(chart_model)
#
#n_patchs = num_cells(chart_model)
#vertex_coords = get_vertex_coordinates(topo)
#TRI_vertex_coords = get_vertex_coordinates(TRI)
#generating_splines = ReferenceFEs._box_splines_222(T, TRI_vertex_coords)
#p_gen_splines = Fill(generating_splines, n_patchs)
#
#tri_to_nodes = get_faces(topo, 2, 0)
#q_vertices = tri_to_nodes[11]
#
#node_adj_graph = compute_graph(topo, 0, 1; self_loops=false)
#vertex_to_vertices = Table(node_adj_graph.rowval, node_adj_graph.colptr)
#
#v_perm = zeros(Int32, 12)
#v_to_v_caches = (array_cache(vertex_to_vertices), array_cache(vertex_to_vertices))
#
#Geometry.compute_loop_patch_vertices!(v_perm, q_vertices, vertex_to_vertices, v_to_v_caches )
#@btime Geometry.compute_loop_patch_vertices!($v_perm, $q_vertices, $vertex_to_vertices, $v_to_v_caches)


# Chart IG fespace
Ω_chart = Triangulation(torus_chart)

T, D = Float64, 3
P = Point{D,Float64}
# Loop reffe has 12 DOF per cell.
loop_reffe = LoopRefFE(P,TRI)
fe = FESpace(Ω_chart, loop_reffe)

alldofs = vcat(fe.cell_dofs_ids...)
@test all( i -> count(==(i), alldofs)==24, alldofs)
@test sort(unique(alldofs)) == 1:(num_vertices(torus_chart)*D)

#Φ = FEFunction(fe, control_vertices)
#err = norm∘(torus_map-Φ)
#writevtk(Ω_chart, "err32", cellfields=["err"=>err]; nsubcells=20)

dΩ = Measure(Ω_chart, 12)
a(Φ_Loop, Ψ_Loop) = ∫(Φ_Loop ⋅ Ψ_Loop)dΩ
l(Ψ_Loop) = ∫(torus_map ⋅ Ψ_Loop)dΩ

op = AffineFEOperator(a, l, fe, fe)

Φh = solve(op)
e = Φh-torus_map
sqrt(sum(a(e,e)))

chart_vertices = get_vertex_coordinates(get_grid_topology(torus_chart))
fitted_control_vertices = evaluate(Φh, chart_vertices )
writevtk(fitted_control_vertices, "fitted_control_vertices")

#freevals = reshape(reinterpret(Float64, hcat(torus_map.(chart_vertices)...)), (600,))
#Φi = FEFunction(fe, freevals)
#recovered_points = evaluate(Φi, chart_vertices) ≠ freevals at all

M = reshape(Φh.free_values, (D, num_vertices(torus_chart)))
fitted_control_vertices = collect( P(r...) for r in eachcol(M))
writevtk(fitted_control_vertices, "fitted_control_vertices")

fitted_loop_torus_model = loop_surface_model(torus_chart, fitted_control_vertices)
writevtk(Triangulation(fitted_loop_torus_model), "fitted_loop_torus"; nsubcells=20)


l2(Ψ_Loop) = ∫(Φh ⋅ Ψ_Loop)dΩ
op2 = AffineFEOperator(a, l2, fe, fe)

Φh2 = solve(op2)
e = Φh-Φh2
@test sqrt(sum(a(e,e))) < 1e-13

end
