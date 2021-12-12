using Gridap.Arrays
using SparseArrays



function shift_to_first(v::Vector, i::T) where {T <: Int}
    circshift(v, -(i - 1))
end

function sort_longest_edge!(elem, node, NT)
    edgelength = zeros(NT, 3)
    #node_v = [[v[1], v[2]] for v in node]
    #node_v = vcat(node_v'...)
    for i = 1:NT
        elem_i = elem[i, :]
        for (j, e) in enumerate(elem_i)
            arr = filter(x -> x != e, elem_i)
            diff = sqrt(sum((node[arr[1]] - node[arr[2]]).^2))
            #diff = sqrt(sum((node_v[arr[1], :] - node_v[arr[2], :]).^2))
            edgelength[i, j] = diff
        end
    end
    max_indices = findmax(edgelength, dims=2)[2]
    for i = 1:NT
        elem[i,:] = shift_to_first(elem[i,:], max_indices[i][2])
    end
    elem
end

function setup_markers_and_nodes!(node, elem, d2p, dualedge, NT, NE, η_arr, θ)
    total = sum(η_arr)
    ix = sortperm(-η_arr)
    current = 0
    marker = zeros(eltype(elem), NE)
    for t = 1:NT
        if (current > θ*total)
            break
        end
        index=1
        ct=ix[t]
        while (index==1)
            base = d2p[elem[ct,2],elem[ct,3]]
            if marker[base]>0
                index=0
            else
                current = current + η_arr[ct]
                N = size(node,1) + 1
                marker[d2p[elem[ct,2],elem[ct,3]]] = N
                midpoint = get_midpoint(node[elem[ct,[2 3],:]])
                node = [node; midpoint]
                ct = dualedge[elem[ct,3],elem[ct,2]]
                if ct==0
                    index=0
                end
            end
        end
    end
    node, marker
end

function divide(elem,t,p)
    elem = [elem; [p[4] p[3] p[1]]]
    elem[t,:]=[p[4] p[1] p[2]]
    elem
end

function bisect(d2p, elem, marker, NT)
    for t=1:NT
        base = d2p[elem[t,2],elem[t,3]]
        if (marker[base]>0)
            p = vcat(elem[t,:], marker[base])
            elem=divide(elem,t,p)
            left=d2p[p[1],p[2]]; right=d2p[p[3],p[1]]
            if (marker[right]>0)
                elem=divide(elem,size(elem,1), [p[4],p[3],p[1],marker[right]])
            end
            if (marker[left]>0)
                elem=divide(elem,t,[p[4],p[1],p[2],marker[left]])
            end
        end
    end
    elem
end

function build_edges(elem)
    edge = [elem[:,[1,2]]; elem[:,[1,3]]; elem[:,[2,3]]]
    unique(sort!(edge, dims=2), dims=1)
end

function build_directed_dualedge(elem, N, NT)
    dualedge = spzeros(eltype(elem), N, N)
    for t=1:NT
        dualedge[elem[t,1],elem[t,2]]=t
        dualedge[elem[t,2],elem[t,3]]=t
        dualedge[elem[t,3],elem[t,1]]=t
    end
    dualedge
end

function dual_to_primal(edge, NE, N)
    d2p = spzeros(eltype(edge), N, N)
    for k=1:NE
        i=edge[k,1]
        j=edge[k,2]
        d2p[i,j]=k
        d2p[j,i]=k
    end
    d2p
end

function test_against_top(face, top::GridTopology, d)
    face_vec = sort.([face[i,:] for i in 1:size(face,1)])
    face_top = sort.(get_faces(top, d, 0))
    issetequal_bitvec = issetequal(face_vec, face_top)
    @assert all(issetequal_bitvec)
end

function get_midpoint(ngon)
    sum(ngon)/length(ngon)
end

function Base.atan(v::VectorValue{2, T}) where {T <: AbstractFloat}
    atan(v[2], v[1])
end

function Base.:^(v::VectorValue{N, T}, r::Core.Integer)  where {T <: AbstractFloat, N}
    VectorValue([v[i]^r for i = 1:N]...)
end

function sort_ccw(cell_coords)
    midpoint = get_midpoint(cell_coords)
    offset_coords = cell_coords .- midpoint
    sorted_perm = sortperm(offset_coords, by=atan)
    return sorted_perm
end

function sort_cell_node_ids_ccw(cell_node_ids, node_coords)
    cell_node_ids_ccw = vcat(cell_node_ids'...)
    #cell_node_ids_ccw = cell_node_ids
    #@show cell_node_ids_ccw
    for (i, cell) in enumerate(cell_node_ids)
        cell_coords = node_coords[cell]
        perm = sort_ccw(cell_coords)
        cell_node_ids_ccw[i,:] = cell[perm]
    end
    cell_node_ids_ccw
end

"""
node_coords == node, cell_node_ids == elem in Long Chen's notation
"""
function newest_vertex_bisection(
    top::GridTopology,
    node_coords::Vector,
    cell_node_ids,
    η_arr,
    θ,
    sort_flag,
)
    # Number of nodes
    N = size(node_coords, 1)
    elem = cell_node_ids
    # Number of cells (triangles in 2D)
    NT = size(elem, 1)
    if sort_flag
        elem = sort_longest_edge!(elem, node_coords, NT)
    end
    # Make sure elem is consistent with GridTopology
    test_against_top(elem, top, 2)
    edge = build_edges(elem)
    NE = size(edge, 1)
    dualedge = build_directed_dualedge(elem, N, NT)
    d2p = dual_to_primal(edge, NE, N)
    # Make sure edge is consistent with GridTopology
    test_against_top(edge, top, 1)
    node_coords, marker =
        setup_markers_and_nodes!(node_coords, elem, d2p, dualedge, NT, NE, η_arr, θ)
    cell_node_ids = bisect(d2p, elem, marker, NT)
    node_coords, cell_node_ids
end

# step 1
function newest_vertex_bisection(
    grid::Grid,
    top::GridTopology,
    η_arr::AbstractVector{<:AbstractFloat},
    θ,
    sort_flag,
)
    node_coords = get_node_coordinates(grid)
    cell_node_ids = get_cell_node_ids(grid)
    if sort_flag
        cell_node_ids_ccw = sort_cell_node_ids_ccw(cell_node_ids, node_coords)
    else
        cell_node_ids_ccw = vcat(cell_node_ids'...)
    end
    node_coords_ref, cell_node_ids_ref =
        newest_vertex_bisection(top, node_coords, cell_node_ids_ccw, η_arr, θ, sort_flag)
    # TODO: Should not convert to matrix and back to Table
    cell_node_ids_ref = Table([c for c in eachrow(cell_node_ids_ref)])
    reffes = get_reffes(grid)
    cell_types = get_cell_type(grid)
    # TODO : Gracefully handle cell_types
    new_cell_types = fill(1, length(cell_node_ids_ref) - length(cell_node_ids))
    append!(cell_types, new_cell_types)
    UnstructuredGrid(node_coords_ref, cell_node_ids_ref, reffes, cell_types)
end

# step 2
function newest_vertex_bisection(
    model::DiscreteModel,
    η_arr::AbstractVector{<:AbstractFloat};
    θ = 1.0, # corresponds to uniform refinement
    sort_flag = false,
)
    grid = get_grid(model)
    top = get_grid_topology(model)
    ref_grid = newest_vertex_bisection(grid, top, η_arr, θ, sort_flag)
    ref_topo = GridTopology(ref_grid)
    #ref_labels = get_face_labeling(model)
    ref_labels = FaceLabeling(ref_topo)
    #ref_labels = # Compute them from the original labels (This is perhaps the most tedious part)
    DiscreteModel(ref_grid, ref_topo, ref_labels)
end