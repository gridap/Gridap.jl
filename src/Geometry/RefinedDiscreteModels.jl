using Gridap.Arrays
using SparseArrays

# TODO: could probably overload ^ from Base for VectorValue
function sumsq(v::VectorValue)
    return sum([v[i]^2 for i in 1:length(v)])
end

function shift_to_first(v::Vector, i::T) where {T <: Int}
    circshift(v, -(i - 1))
end

# TODO: under construction for longest side in non-uniform mesh
function sort_longest_side(node::Vector, elem::Matrix, NT)
    edgelength = zeros(NT, 3)
    node = [[v[1], v[2]] for v in node]
    node = vcat(node'...)
    for i = 1:NT
        elem_i = elem[i, :]
        for (j, e) in enumerate(elem_i)
            arr = filter(x -> x != e, elem_i)
            diff = sqrt(sum((node[arr[1], :] - node[arr[2], :]).^2))
            edgelength[i, j] = diff
        end
    end
    max_indices = findmax(edgelength, dims=2)[2]
    for i = 1:NT
        max_indices[i][2]
        elem[i,:] = shift_to_first(elem[i,:], max_indices[i][2])
    end
    elem
end

function setup_markers(NT, NE, node, elem, d2p, dualedge, θ)
    # TODO: handle estimator
    η = fill(1, NT)
    total = sum(η)
    ix = sortperm(-η)
    current = 0
    marker = zeros(Int32, NE,1)
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
                current = current + η[ct]
                N = size(node,1)+1;
                marker[d2p[elem[ct,2],elem[ct,3]]] = N
                midpoint = get_midpoint(node[elem[ct,[2 3],:]])
                node = [node; midpoint]
                ct = dualedge[elem[ct,2],elem[ct,3]]
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
    # TODO: for now everything marked for bisection
    for t=1:NT
        base=d2p[elem[t,2],elem[t,3]]
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

function build_edges(elem::Matrix{T}) where {T <: Integer}
    edge = [elem[:,[1,2]]; elem[:,[1,3]]; elem[:,[2,3]]]
    unique(sort!(edge, dims=2), dims=1)
end

function build_directed_dualedge(elem::Matrix{Ti}, N::T, NT::T) where {Ti, T <: Integer}
    dualedge = spzeros(Int64, N, N)
    for t=1:NT
        dualedge[elem[t,1],elem[t,2]]=t
        dualedge[elem[t,2],elem[t,3]]=t
        dualedge[elem[t,3],elem[t,1]]=t
    end
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
    face_vec = sort.([face[i,:] for i in 1:size(face,1)])
    face_top = sort.(get_faces(top, d, 0))
    issetequal_bitvec = issetequal(face_vec, face_top)
    @assert all(issetequal_bitvec)
end

function get_midpoint(ngon)
    return sum(ngon)/length(ngon)
end

function Base.atan(v::VectorValue{2, T}) where {T <: AbstractFloat}
    atan(v[2], v[1])
end

function sort_ccw(cell_coords)
    midpoint = get_midpoint(cell_coords)
    offset_coords = cell_coords .- midpoint
    sorted_perm = sortperm(offset_coords, by=atan)
    return sorted_perm
end

function sort_cell_node_ids_ccw(cell_node_ids, node_coords)
    cell_node_ids_ccw = vcat(cell_node_ids'...)
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
function newest_vertex_bisection(top::GridTopology, node_coords::Vector, cell_node_ids::Matrix)
    N = size(node_coords, 1)
    #@show elem = vcat(cell_node_ids'...)
    elem = cell_node_ids
    NT = size(elem, 1)
    elem = sort_longest_side(node_coords, elem, NT)
    test_against_top(elem, top, 2)
    edge = build_edges(elem)
    NE = size(edge, 1)
    dualedge = build_directed_dualedge(elem, N, NT)
    d2p = dual_to_primal(edge, NE, N)
    test_against_top(edge, top, 1)
    ## TODO: Mark largest edge
    ##sort_elem_for_labeling(node_coords, elem)
    node_coords, marker = setup_markers(NT, NE, node_coords, elem, d2p, dualedge, 1)
    #@show node
    @show size(node_coords, 1)
    cell_node_ids_ref = bisect(d2p, elem, marker, NT)
    @show size(cell_node_ids_ref, 1)
    #@show cell_node_ids_ref
    node_coords, cell_node_ids_ref
end

# step 1
function newest_vertex_bisection(grid::Grid, top::GridTopology, cell_mask::AbstractVector{<:Bool})
    #get_faces(top
    #@show cell_coords = get_cell_coordinates(grid)
    node_coords = get_node_coordinates(grid)
    cell_node_ids = get_cell_node_ids(grid)
    cell_node_ids_ccw = sort_cell_node_ids_ccw(cell_node_ids, node_coords)
    typeof(node_coords)
    # TODO: Modify node__coords and cell_node_ids
    #node_coords, cell_node_ids = newest_vertex_bisection(top, node_coords, cell_node_ids_ccw)
    newest_vertex_bisection(top, node_coords, cell_node_ids_ccw)
    reffes = get_reffes(grid)
    cell_types = get_cell_type(grid)
    UnstructuredGrid(node_coords, cell_node_ids, reffes, cell_types)
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
