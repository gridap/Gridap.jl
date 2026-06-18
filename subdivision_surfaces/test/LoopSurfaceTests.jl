module LoopSurfaceTests

using Test
using Gridap
using GridapSubdivisionSurfaces
using GridapSubdivisionSurfaces.ReferenceFEs

using Gridap.ReferenceFEs
using Gridap.Fields

using Gridap.Geometry
using Gridap.Geometry: get_patch_faces
using Gridap.Geometry: get_patch_cells
using Gridap.Geometry: get_patch_to_tfaces
using Gridap.Geometry: get_pface_to_lpface
using Gridap.Geometry: get_pface_to_patch
using Gridap.Geometry: num_patches
using Gridap.Geometry: extend_patches_by_single_layer

"""
    loop_torus_chart_and_map(n1=5, n2=2*n1)

Torus chart grid ([0,1]²) meshed with 2*n1*n2 triangles.
"""
function loop_torus_chart_and_map(n1=5, n2=2*n1)
  model = CartesianDiscreteModel((0,1,0,1), (n1,n2), isperiodic=(true,true))
  model = simplexify(model)
  torus_chart = Gridap.Geometry.UnstructuredDiscreteModel(model)

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

torus_chart, torus_map = loop_torus_chart_and_map(5)
loop_torus_model = loop_surface_model(torus_chart, torus_map)
loop_torus_trian = Triangulation(loop_torus_model)
writevtk(loop_torus_trian, "loop_torus"; nsubcells=20)

control_vertices = get_vertex_coordinates(get_grid_topology(loop_torus_model ));
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

alldofs = vcat(fe.cell_dofs_ids...);
@test all( i -> count(==(i), alldofs)==24, alldofs)
@test sort(unique(alldofs)) == 1:(num_vertices(torus_chart)*D)

dΩ = Measure(Ω_chart, 8)
a(Φ_Loop, Ψ_Loop) = ∫(Φ_Loop ⋅ Ψ_Loop)dΩ
l(Ψ_Loop) = ∫(torus_map ⋅ Ψ_Loop)dΩ

op = AffineFEOperator(a, l, fe, fe)

Φh = solve(op)
e = Φh-torus_map
err = sqrt(sum(a(e,e)))

#chart_vertices = get_vertex_coordinates(get_grid_topology(torus_chart));
#fitted_torus_vertices = evaluate(Φh, chart_vertices);
#writevtk(fitted_torus_vertices, "fitted_torus_vertices")

M = reshape(Φh.free_values, (D, num_vertices(torus_chart)));
fitted_control_vertices = collect( P(r...) for r in eachcol(M));
writevtk(fitted_control_vertices, "fitted_control_vertices")

topo = get_grid_topology(torus_chart)
fitted_loop_torus_model = loop_surface_model(topo , fitted_control_vertices)
writevtk(Triangulation(fitted_loop_torus_model), "fitted_loop_torus"; nsubcells=20)

# H1 projection
chart_model = torus_chart
geo_map = torus_map
# get (loosely) approximate control vertex coordinates by interpolating the geo_map
geo_map_field = GenericField(geo_map)
chart_vertex_coordinates = get_vertex_coordinates(get_grid_topology(chart_model))
geo_map_cache = return_cache(geo_map_field, chart_vertex_coordinates)
vertex_coordinates = evaluate!(geo_map_cache, geo_map_field, chart_vertex_coordinates)

# Definition of a Loop FE space
P = eltype(vertex_coordinates)
loop_reffe = LoopRefFE(P,TRI)
Ω_chart = Triangulation(chart_model)
dΩ = Measure(Ω_chart, 8)
dΩf = Measure(Ω_chart, 16)
l2(e) = sqrt(sum( ∫( e⋅e)dΩf ))
h1(e) = sqrt(sum( ∫(∇(e)⊙∇(e) + e⋅e)dΩf ))
# this dispatches to hacky subdivision_surfaces/src/FESpaces/LoopFESpaces.jl
U = FESpace(Ω_chart, loop_reffe)

# L2 projection
a_l2(Φ, Ψ) = ∫(Ψ⋅Φ)dΩ
l(Ψ) = ∫(Ψ⋅geo_map)dΩ
op_l2 = AffineFEOperator(a_l2, l, U, U)
Φh_l2 = solve(op_l2)
e_l2 = Φh_l2-geo_map

l2P_errl2 = l2(e_l2)
l2P_errh1 = h1(e_l2)


# H1 projection
ε = 1.e-0
#ε = 1.e-6
a_h1(Φ, Ψ) = ∫(∇(Ψ)⊙∇(Φ) + ε*Ψ⋅Φ)dΩ
l(   Ψ) = ∫(∇(Ψ)⊙∇(geo_map) + ε*Ψ⋅geo_map)dΩ
op_h1 = AffineFEOperator(a_h1, l, U, U)
Φh_h1 = solve(op_h1)
e_h1 = Φh_h1-geo_map

h1P_errl2 = l2(e_h1)
h1P_errh1 = h1(e_h1)

@info "Loop surface L²-projection error: L² error: $l2P_errl2, H¹ error: $l2P_errh1 "
@info "Loop surface H¹-projection error: L² error: $h1P_errl2, H¹ error: $h1P_errh1"

# Retrieving the control vertices coordinates from the DoF vector
D = length(P)
n_vertices = length(vertex_coordinates)
M = reshape(Φh.free_values, (D, n_vertices));
fitted_control_vertices = collect( P(r...) for r in eachcol(M));

loop_surface_model(get_grid_topology(chart_model), fitted_control_vertices)

end
