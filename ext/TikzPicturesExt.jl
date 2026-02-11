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
- `label_offset`: offset for labels, w.r.t. nodes (default: 0.2)
- `edge_style`: TikZ style for edges (default: "[ ]")
- `node_style`: TikZ style for nodes (default: "[circle, fill=black, inner sep=0pt, minimum size=0.1cm]")
- `label_style`: TikZ style for labels (default: "[ ]")

as well as any `TikzPicture`-specific keyword arguments.

"""
function TikzPictures.TikzPicture(model::DiscreteModel;kwargs...)
  TikzPicture(get_grid(model);kwargs...)
end

function TikzPictures.TikzPicture(grid::Grid;kwargs...)
  TikzPicture(GridTopology(grid);kwargs...)
end

function TikzPictures.TikzPicture(topo::GridTopology;kwargs...)
  edge_to_nodes = Geometry.get_faces(topo, 1, 0)
  node_coordinates = Geometry.get_vertex_coordinates(topo)
  draw_graph(edge_to_nodes, node_coordinates; kwargs...)
end

function TikzPictures.TikzPicture(poly::Polytope;kwargs...)
  edge_to_nodes = Geometry.get_faces(poly, 1, 0)
  node_coordinates = Geometry.get_vertex_coordinates(poly)
  draw_graph(edge_to_nodes, node_coordinates; kwargs...)
end

function draw_graph(
  edge_to_nodes, node_coordinates;
  draw_nodes = true,
  draw_labels = true,
  label_offset = 0.2,
  edge_style = "[ ]",
  node_style = "[circle, fill=black, inner sep=0pt, minimum size=0.1cm]",
  label_style = "[ ]",
  tikz_kwargs...
)
  n_edges = length(edge_to_nodes)
  n_nodes = length(node_coordinates)

  draw_pt(v) = "(" * join(Tuple(v),",") * ")"
  draw_edge(e) = draw_edge(edge_to_nodes[e]...)
  draw_edge(a,b) = "\\draw $(edge_style)" * draw_pt(node_coordinates[a]) * " -- " * draw_pt(node_coordinates[b]) * "; \n"
  draw_node(a) = "\\node $(node_style) at " * draw_pt(node_coordinates[a]) * " {}; \n"
  draw_label(a) = "\\node $(label_style) at " * draw_pt(node_coordinates[a] .+ label_offset) * " { $a }; \n"

  s = join((draw_edge(e) for e in 1:n_edges), "")
  if draw_nodes
  s *= join((draw_node(a) for a in 1:n_nodes), "")
  end
  if draw_labels
    s *= join((draw_label(a) for a in 1:n_nodes), "")
  end
  return TikzPicture(s;tikz_kwargs...)
end

end