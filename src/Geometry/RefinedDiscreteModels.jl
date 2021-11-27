using Gridap.Arrays

function get_midpoint(x::AbstractVector, y::AbstractVector)
    (x + y) ./ 2.0
end

function build_edges(elem::Matrix{T}) where {T <: Integer}
    edge = [elem[:,[1,2]]; elem[:,[1,3]]; elem[:,[2,3]]]
    unique(sort!(edge, dims=2), dims=1)
end

function build_dual_edge(elem::Matrix{Ti}, N::T) where {Ti, T <: Integer}
    NT = size(elem, 1)
    [elem[:,[1,2,3]],elem[:,[2,3,1]], [1:NT,1:NT,1:NT],N,N]
end

function test_against_top(edge::Matrix{Ti}, top::GridTopology, d::T) where {Ti, T <: Integer}
    edge_vec = [edge[i,:] for i in 1:size(edge,1)]
    edge_top = get_faces(top, d, 0)
    @assert issetequal(edge_vec, edge_top)
end

"""
node_coords == node, cell_node_ids == elem in Long Chen's notation
"""
function newest_vertex_bisection(top::GridTopology, node_coords::Vector, cell_node_ids::Table)
    N = size(node_coords, 1)
    elem = vcat(cell_node_ids'...)
    test_against_top(elem, top, 2)
    edge = build_edges(elem)
    @show dual_edge = build_dual_edge(elem, N)
    test_against_top(edge, top, 1)
    node_coords, cell_node_ids
end

# setp 1
function newest_vertex_bisection(grid::Grid, top::GridTopology, cell_mask::AbstractVector{<:Bool})
    #get_faces(top
    #@show cell_coords = get_cell_coordinates(grid)
    node_coords = get_node_coordinates(grid)
    cell_node_ids = get_cell_node_ids(grid)
    # TODO: Modify node__coords and cell_node_ids
    typeof(node_coords)
    node_coords, cell_node_ids = newest_vertex_bisection(top, node_coords, cell_node_ids)
    reffes = get_reffes(grid)
    cell_types = get_cell_type(grid)
    UnstructuredGrid(node_coords, cell_node_ids, reffes, cell_types)
    #@show new_val = VectorValue{2, Float64}(1.5, 1.5, 1.5)
    #new_val = LazyArray([VectorValue(1.5, 1.5), VectorValue(1.2, 2), VectorValue(3.2, 3.1)])
    #@show new_val = lazy_map(get_midpoint, cell_coords[1:2], 0*cell_coords[1:2])
    #@show new_val = lazy_map(get_midpoint, cell_coords[1:2], 0*cell_coords[1:2])
    #cell_coords = lazy_append(cell_coords, new_val)
  # tod
  #ref_grid
end

# step 2
function newest_vertex_bisection(model::DiscreteModel,cell_mask::AbstractVector{<:Bool})
  grid  = get_grid(model)
  top = get_grid_topology(model)
  ref_grid = newest_vertex_bisection(grid, top, cell_mask)
  #ref_topo = GridTopology(grid)
  #labels = get_face_labelling(model)
  #ref_labels = # Compute them from the original labels (This is perhaps the most tedious part)
  #ref_model = DiscreteModel(ref_grid,ref_topo,ref_labels)
  #ref_model
end
