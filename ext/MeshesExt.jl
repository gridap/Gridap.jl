module MeshesExt

using Meshes: Meshes
using Gridap

using Gridap.Geometry
using Gridap.Geometry: Grid, get_node_coordinates, get_cell_node_ids
using Gridap.Geometry: GridTopology, get_vertex_coordinates, get_faces

function Meshes.viz(grid::Grid;kwargs...)
  nodes = get_node_coordinates(grid)
  connectivity = get_cell_node_ids(grid)
  _viz(nodes,connectivity;kwargs...)
end

function Meshes.viz(topo::GridTopology{Dc};kwargs...) where Dc
  nodes = get_vertex_coordinates(topo)
  connectivity = get_faces(topo,Dc,0)
  _viz(nodes,connectivity;kwargs...)
end

function _viz(nodes,connectivity;kwargs...)
  perms = map(n -> orient_nodes(nodes[n]), connectivity)
  elems = map((n,p) -> Meshes.connect(Tuple(n[p])), connectivity, perms)
  verts = map(v -> Meshes.Point(v...), nodes)
  mesh = Meshes.SimpleMesh(verts, elems)
  Meshes.viz(mesh;kwargs...)
end

function orient_nodes(v)
  _angle(c,v) = atan(v[1]-c[1],v[2]-c[2])
  c = mean(v)
  sortperm(v, by = x -> _angle(c,x), rev = true)
end

end