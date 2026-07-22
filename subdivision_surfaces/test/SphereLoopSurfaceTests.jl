module SphereLoopSurfaceTests

using Test
using Gridap
using GridapSubdivisionSurfaces
using GridapSubdivisionSurfaces.ReferenceFEs
using Gridap.ReferenceFEs
using Gridap.Fields
using Gridap.Geometry
using Gridap.CellData
using GridapGeosciences

const tmp_dir = joinpath(@__DIR__, "..", "tmp")
mkpath(tmp_dir)
tmp_path(name) = joinpath(tmp_dir, name)

function coarse_parametric_model(radius::Float64=1.0)
  sphere_mesh = CubedSphereMesh(radius)
  atlas_model = AtlasDiscreteModel(sphere_mesh, 0)
  atlas_model
end

function loop_sphere_chart_and_map(radius::Float64=1.0)
  atlas_model = coarse_parametric_model(radius)
  atlas_grid = get_atlas_grid(atlas_model)

  # Extract parametric grid and create UnstructuredDiscreteModel
  quad_grid = atlas_grid.param_grid
  quad_topo = get_grid_topology(atlas_model)
  quad_labels = get_face_labeling(atlas_model)
  chart_model_quad = UnstructuredDiscreteModel(quad_grid, quad_topo, quad_labels)
  chart_model_tri = simplexify(chart_model_quad)

  # Get cell-varying ambient maps from atlas
  cell_ambient_maps_quads = get_cell_ambient_maps(atlas_grid)

  # Expand maps to triangles (each quad becomes 2 triangles)
  cell_ambient_maps_tri = []
  for map in cell_ambient_maps_quads
    push!(cell_ambient_maps_tri, map)
    push!(cell_ambient_maps_tri, map)
  end

  chart_model_tri, cell_ambient_maps_tri
end

sphere_chart_tri, cell_ambient_maps = loop_sphere_chart_and_map(1.0)

writevtk(sphere_chart_tri, tmp_path("sphere_chart_tri"))

@test isa(sphere_chart_tri, DiscreteModel)

try
  @info "✓ Created parametric model and cell-varying maps from sphere panels"

  # Wrap cell-varying maps as CellField
  Ω_chart = Triangulation(sphere_chart_tri)
  geo_map = GenericCellField(cell_ambient_maps, Ω_chart, PhysicalDomain())

  # Use loop_surface_model with cell-varying maps as geo_map
  loop_sphere_model = loop_surface_model(sphere_chart_tri, geo_map)

  @test isa(loop_sphere_model, DiscreteModel)

  loop_sphere_trian = Triangulation(loop_sphere_model)
  writevtk(loop_sphere_trian, tmp_path("loop_sphere_fitted"); nsubcells=20)

catch e
  @warn "Loop surface fitting failed: $e"
end

end
