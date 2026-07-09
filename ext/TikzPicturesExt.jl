module TikzPicturesExt

using TikzPictures
using Gridap
using Gridap.Geometry

"""
    TikzPictures.TikzPicture(model::DiscreteModel;kwargs...)
    TikzPictures.TikzPicture(grid::Grid;kwargs...)
    TikzPictures.TikzPicture(topo::GridTopology;kwargs...)
    TikzPictures.TikzPicture(topo::Polytope;kwargs...)

Returns a `TikzPicture` for the given model, grid, or topology. 

Available keyword arguments include:

- `draw_nodes`: whether to draw visible nodes (default: true)
- `draw_labels`: whether to add labels (default: true)
- `edge_style`: TikZ style for edges (default: "[ ]")
- `node_style`: TikZ style for nodes (default: "[circle, fill=black, inner sep=0pt, minimum size=0.1cm]")
- `label_style`: TikZ style for labels (default: "[ ]")
- `label_offset`: offset for labels, w.r.t. nodes (default: 0.0)
- `label_dim`: dimension of the labeled faces (default: 0, i.e nodes)

as well as any `TikzPicture`-specific keyword arguments.

"""
function TikzPictures.TikzPicture(model::DiscreteModel;kwargs...)
  TikzPicture(get_grid(model);kwargs...)
end

function TikzPictures.TikzPicture(grid::Grid;kwargs...)
  TikzPicture(GridTopology(grid);kwargs...)
end

function TikzPictures.TikzPicture(topo::GridTopology;kwargs...)
  draw_topology(topo; kwargs...)
end

function TikzPictures.TikzPicture(poly::Polytope;kwargs...)
  draw_topology(poly; kwargs...)
end

function TikzPictures.TikzPicture(rr::Gridap.Adaptivity.RefinementRule;kwargs...)
  TikzPicture(rr.reg_grid; kwargs...)
end

function draw_topology(
  topo;
  draw_nodes = true,
  draw_labels = true,
  edge_style = "[ ]",
  node_style = "[circle, fill=black, inner sep=0pt, minimum size=0.1cm]",
  label_style = "[ ]",
  label_dim = 0,
  label_offset = 0.0,
  tikz_kwargs...
)
  n_nodes = num_vertices(topo)
  n_edges = num_faces(topo, 1)
  edge_to_nodes = Geometry.get_faces(topo, 1, 0)
  node_coordinates = Geometry.get_vertex_coordinates(topo)

  draw_pt(v) = "(" * join(Tuple(v),",") * ")"
  draw_edge(e) = draw_edge(edge_to_nodes[e]...)
  draw_edge(a,b) = "\\draw $(edge_style)" * draw_pt(node_coordinates[a]) * " -- " * draw_pt(node_coordinates[b]) * "; \n"
  draw_node(a) = "\\node $(node_style) at " * draw_pt(node_coordinates[a]) * " {}; \n"
  draw_label(f,v) = "\\node $(label_style) at " * draw_pt(v) * " { $f }; \n"

  s = join((draw_edge(e) for e in 1:n_edges), "")
  if draw_nodes
    s *= join((draw_node(a) for a in 1:n_nodes), "")
  end
  if draw_labels
    face_to_nodes = Geometry.get_faces(topo, label_dim, 0)
    face_centroids = map(nodes -> mean(node_coordinates[nodes]) .+ label_offset, face_to_nodes)
    s *= join((draw_label(f,v) for (f,v) in enumerate(face_centroids)), "")
  end
  return TikzPicture(s;tikz_kwargs...)
end

end