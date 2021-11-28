using Gridap.Arrays
using SparseArrays
using Printf

function build_edges(elem::Matrix{T}) where {T <: Integer}
    edge = [elem[:,[1,2]]; elem[:,[1,3]]; elem[:,[2,3]]]
    unique(sort!(edge, dims=2), dims=1)
end

function build_directed_dual_edge(elem::Matrix{Ti}, N::T) where {Ti, T <: Integer}
    NT = size(elem, 1)
    # TODO: Make sparse?
    #@show elem[:,[1,2,3]]
    #sparse(elem[:,[1,2,3]],elem[:,[2,3,1]], [1:NT,1:NT,1:NT])
    #elem[2, :] = [2, 4, 3]
    dualedge = spzeros(Int64, N, N)
    for t=1:NT
        #@show t
        #@show elem[t,1],elem[t,2]
        dualedge[elem[t,1],elem[t,2]]=t
        #show elem[t,2],elem[t,3]
        dualedge[elem[t,2],elem[t,3]]=t
        #@show elem[t,3],elem[t,1]
        dualedge[elem[t,3],elem[t,1]]=t
    end
    #@show nnz(dualedge)
    dualedge
end

function dual_to_primal(edge::Matrix{Ti}, NE::T, N::T) where {Ti, T <: Integer}
    d2p = spzeros(Int64, N, N)
    for k=1:NE
        i=edge[k,1]
        j=edge[k,2]
        d2p[i,j]=k
        d2p[j,i]=k
    end
    d2p
end

function test_against_top(face::Matrix{Ti}, top::GridTopology, d::T) where {Ti, T <: Integer}
    face_vec = [face[i,:] for i in 1:size(face,1)]
    @show face_top = get_faces(top, d, 0)
    @show face_top
    @assert issetequal(face_vec, face_top)
end

"""
node_coords == node, cell_node_ids == elem in Long Chen's notation
"""
function newest_vertex_bisection(top::GridTopology, node_coords::Vector, cell_node_ids::Table)
    N = size(node_coords, 1)
    @show elem = vcat(cell_node_ids'...)
    test_against_top(elem, top, 2)
    edge = build_edges(elem)
    NE = size(edge, 1)
    @show dual_edge = build_directed_dual_edge(elem, N)
    d2p = dual_to_primal(edge, NE, N)
    test_against_top(edge, top, 1)
    node_coords, cell_node_ids
end


function get_midpoint(ngon::Vector)
    return sum(ngon)/length(ngon)
end

function Base.atan(v::VectorValue{2, T}) where {T <: AbstractFloat}
    atan(v[2], v[1])
end

function sort_ccw(cell_coords)
    midpoint = get_midpoint(cell_coords)
    offset_coords = cell_coords .- midpoint
    #@show midpoint_zero = get_midpoint(offset_coords)
    sorted_perm = sortperm(offset_coords, by=atan)
    @show offset_coords[sorted_perm] .+ midpoint
    return sorted_perm
end


function resort_cell_node_ids(cell_node_ids, node_coords)
    for cell in cell_node_ids
        @show cell_coords = node_coords[cell]
        @show perm = sort_ccw(cell_coords)
        @show cell[perm]
    end
end

# setp 1
function newest_vertex_bisection(grid::Grid, top::GridTopology, cell_mask::AbstractVector{<:Bool})
    #get_faces(top
    #@show cell_coords = get_cell_coordinates(grid)
    node_coords = get_node_coordinates(grid)
    cell_node_ids = get_cell_node_ids(grid)
    ccw_cell_node_ids = resort_cell_node_ids(cell_node_ids, node_coords)
    # TODO: Modify node__coords and cell_node_ids
    typeof(node_coords)
    #node_coords, cell_node_ids = newest_vertex_bisection(top, node_coords, cell_node_ids)
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
