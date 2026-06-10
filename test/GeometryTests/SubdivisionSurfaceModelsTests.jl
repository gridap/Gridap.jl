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
function simplex_torus_chart(n=4)
  model = CartesianDiscreteModel((0,1,0,1), (n,n), isperiodic=(true,true))
  model = simplexify(model)
  UnstructuredDiscreteModel(model)
end
torus_chart = simplex_torus_chart(10)
writevtk(torus_chart , "torus_chart")

function torus_map(p)
  x,y = p
  R = 5
  r = 1
  θ = 2π*x
  φ = 2π*y
  Point( (R+r*sin(θ))*cos(φ), (R+r*sin(θ))*sin(φ), r*cos(θ) )
end

loop_torus_model = loop_surface_model(torus_chart, torus_map)
loop_torus_trian = Triangulation(loop_torus_model)
writevtk(loop_torus_trian, "loop_torus_model"; nsubcells=20)
