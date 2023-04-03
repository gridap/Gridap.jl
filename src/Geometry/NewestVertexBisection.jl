# This implementation is based on the following article:
#
# LONG CHEN (2008). SHORT IMPLEMENTATION OF BISECTION IN MATLAB.
# In Recent Advances in Computational Sciences. WORLD SCIENTIFIC.
# DOI: 10.1142/9789812792389_0020


using SparseArrays
using Random
using AbstractTrees

# The following was taken directly from Tim Holy's example at
# https://github.com/JuliaCollections/AbstractTrees.jl/blob/master/examples/binarytree_core.jl
#####################################################################################
mutable struct BinaryNode{T}
  data::T
  parent::BinaryNode{T}
  left::BinaryNode{T}
  right::BinaryNode{T}

  # Root constructor
  BinaryNode{T}(data) where {T} = new{T}(data)
  # Child node constructor
  BinaryNode{T}(data, parent::BinaryNode{T}) where {T} = new{T}(data, parent)
end
BinaryNode(data) = BinaryNode{typeof(data)}(data)

function leftchild(data, parent::BinaryNode)
  !isdefined(parent, :left) || error("left child is already assigned")
  node = typeof(parent)(data, parent)
  parent.left = node
end
function rightchild(data, parent::BinaryNode)
  !isdefined(parent, :right) || error("right child is already assigned")
  node = typeof(parent)(data, parent)
  parent.right = node
end

function AbstractTrees.children(node::BinaryNode)
  if isdefined(node, :left)
    if isdefined(node, :right)
      return (node.left, node.right)
    end
    return (node.left,)
  end
  isdefined(node, :right) && return (node.right,)
  return ()
end

AbstractTrees.printnode(io::IO, node::BinaryNode) = print(io, node.data)

Base.eltype(::Type{<:TreeIterator{BinaryNode{T}}}) where {T} = BinaryNode{T}
Base.IteratorEltype(::Type{<:TreeIterator{BinaryNode{T}}}) where {T} = Base.HasEltype()
#####################################################################################

_are_parallel(v, w) = v[1] * w[2] == v[2] * w[1]

_get_midpoint(ngon::AbstractArray{<:VectorValue}) = sum(ngon) / length(ngon)

_shift_to_first(v::AbstractVector{T}, i::T) where {T<:Integer} = circshift(v, -(i - 1))

function _print_forest(forest::AbstractArray{<:BinaryNode})
  num_leaves = 0
  println("_____________________")
  for root_cell in forest
    print_tree(root_cell)
    num_leaves += length(collect(Leaves(root_cell)))
  end
  println("_____________________")
  @show num_leaves
end

function _set_d_to_dface_to_old_node!(
  d_to_dface_to_oldid::AbstractArray{<:AbstractArray},
  d_to_dface_to_olddim::AbstractArray{<:AbstractArray},
  topo::GridTopology,
  node_to_bis_edge::Dict,
  num_new_nodes::Integer,
)

  edge_to_node = get_faces(topo, 1, 0)
  push!(d_to_dface_to_oldid, Vector{}(undef, num_new_nodes))
  push!(d_to_dface_to_olddim, Vector{}(undef, num_new_nodes))
  # Loop over all new node IDs
  for node_id = 1:num_new_nodes
    # New node bisected an edge
    if node_id in keys(node_to_bis_edge)
      bis_edge = node_to_bis_edge[node_id]
      bis_edge_index = findfirst(edge -> edge == bis_edge, edge_to_node)
      d_to_dface_to_oldid[1][node_id] = bis_edge_index
      d_to_dface_to_olddim[1][node_id] = 1
    # Node corresponds to previous node with same ID
    # since node IDs are appended to the end
    else
      d_to_dface_to_oldid[1][node_id] = node_id
      d_to_dface_to_olddim[1][node_id] = 0
    end
  end
end

function _set_d_to_dface_to_old_edge!(
  d_to_dface_to_oldid::AbstractArray{<:AbstractArray},
  d_to_dface_to_olddim::AbstractArray{<:AbstractArray},
  forest::AbstractArray{<:BinaryNode},
  topo::GridTopology,
  topo_ref::GridTopology,
  vertices::AbstractArray{<:VectorValue},
)
  # Topological mappings for the refined mesh
  edge_to_node_ref = get_faces(topo_ref, 1, 0)
  cell_to_edge_ref = get_faces(topo_ref, 2, 1)
  num_edges = length(edge_to_node_ref)
  # Topological mappings for the old mesh
  cell_to_edge = get_faces(topo, 2, 1)
  edge_to_node = get_faces(topo, 1, 0)
  cell_to_node = get_faces(topo, 2, 0)
  push!(d_to_dface_to_oldid, Vector{}(undef, num_edges))
  push!(d_to_dface_to_olddim, Vector{}(undef, num_edges))
  for root_cell in forest
    root_edges = cell_to_edge[root_cell.data]
    root_nodes = cell_to_node[root_cell.data]
    for leaf_cell in Leaves(root_cell)
      # Set because of possible repeats for new cells that are neighbors
      leaf_edges = Set(cell_to_edge_ref[leaf_cell.data])
      for leaf_edge in leaf_edges
        leaf_edge_nodes = edge_to_node_ref[leaf_edge]
        is_on_top_of_former_edge = false
        # We need to check if this edge is parallel to any of the
        # other edges in the old cell to see if it is on top of a former
        # edge
        if !isempty(intersect(leaf_edge_nodes, root_nodes))
          for root_edge in root_edges
            root_edge_nodes = edge_to_node[root_edge]
            root_edge_vector = vertices[root_edge_nodes[2]] - vertices[root_edge_nodes[1]]
            leaf_edge_vector = vertices[leaf_edge_nodes[2]] - vertices[leaf_edge_nodes[1]]
            if _are_parallel(root_edge_vector, leaf_edge_vector)
              d_to_dface_to_oldid[2][leaf_edge] = root_edge
              d_to_dface_to_olddim[2][leaf_edge] = 1
              is_on_top_of_former_edge = true
            end
          end
        end
        # Default case: this edge bisects a cell, i.e. not on top of
        # former edge
        if !is_on_top_of_former_edge
          d_to_dface_to_oldid[2][leaf_edge] = root_cell.data
          d_to_dface_to_olddim[2][leaf_edge] = 2
        end
      end
    end
  end
end

function _set_d_to_dface_to_old_cell!(
  d_to_dface_to_oldid::AbstractArray{<:AbstractArray},
  d_to_dface_to_olddim::AbstractArray{<:AbstractArray},
  forest::AbstractArray{<:BinaryNode},
  topo_ref::GridTopology,
)
  cell_to_node_ref = get_faces(topo_ref, 2, 0)
  num_cells = length(cell_to_node_ref)
  push!(d_to_dface_to_oldid, Vector{}(undef, num_cells))
  push!(d_to_dface_to_olddim, Vector{}(undef, num_cells))
  # Can basically just use the forest. All new cells should point
  # to the root of their tree
  leaf_ids = Int32[]
  for root_cell in forest
    root_id = root_cell.data
    for leaf_cell in Leaves(root_cell)
      leaf_id = leaf_cell.data
      push!(leaf_ids, leaf_id)
      d_to_dface_to_oldid[3][leaf_id] = root_id
      d_to_dface_to_olddim[3][leaf_id] = 2
    end
  end
end

function _create_d_to_dface_to_old(
  forest::AbstractArray{<:BinaryNode},
  topo::GridTopology,
  topo_ref::GridTopology,
  vertices::AbstractArray{<:VectorValue},
  node_to_bis_edge::Dict,
)
  d_to_dface_to_oldid = Vector{Vector}()
  d_to_dface_to_olddim = Vector{Vector}()
  num_new_nodes = length(vertices)
  # The order is important here:
  # Node
  _set_d_to_dface_to_old_node!(
    d_to_dface_to_oldid,
    d_to_dface_to_olddim,
    topo,
    node_to_bis_edge,
    num_new_nodes,
  )
  # Edge
  _set_d_to_dface_to_old_edge!(
    d_to_dface_to_oldid,
    d_to_dface_to_olddim,
    forest,
    topo,
    topo_ref,
    vertices,
  )
  # Cell
  _set_d_to_dface_to_old_cell!(
     d_to_dface_to_oldid,
     d_to_dface_to_olddim,
     forest, topo_ref)
  # Checks to make sure all new dims and ids are filled
  for d = 0:2
    @test undef ∉ d_to_dface_to_oldid[d+1]
    @test undef ∉ d_to_dface_to_olddim[d+1]
    @test length(d_to_dface_to_olddim[d+1]) == length(d_to_dface_to_olddim[d+1])
  end
  d_to_dface_to_olddim, d_to_dface_to_oldid
end

function _propagate_labeling(model, d_to_dface_to_olddim, d_to_dface_to_oldid)
  labels = get_face_labeling(model)
  labels_ref = FaceLabeling(length.(d_to_dface_to_oldid))
  @test isempty(labels_ref.tag_to_entities)
  # Copy entities
  for entities in labels.tag_to_entities
    push!(labels_ref.tag_to_entities, copy(entities))
  end
  # Copy refs
  @test isempty(labels_ref.tag_to_name)
  for name in labels.tag_to_name
    push!(labels_ref.tag_to_name, name)
  end
  # Populate new entity IDs using d_to_dface_to_oldid and
  # d_to_dface_to_olddim
  for d = 0:num_dims(model)
    oldid_d = d_to_dface_to_oldid[d+1]
    olddim_d = d_to_dface_to_olddim[d+1]
    @test length(oldid_d) == length(olddim_d)
    for i = 1:length(oldid_d)
      old_id = oldid_d[i]
      old_dim = olddim_d[i]
      old_entity_id = labels.d_to_dface_to_entity[old_dim+1][old_id]
      labels_ref.d_to_dface_to_entity[d+1][i] = old_entity_id
    end
  end
  labels_ref
end

function _sort_longest_edge!(
  elem::AbstractVector{<:AbstractVector{T}},
  node::AbstractVector{<:VectorValue},
) where {T<:Integer}
  NT = length(elem)
  edgelength = zeros(NT, 3)
  for t = 1:NT
    elem_t = elem[t][:]
    for (j, e) in enumerate(elem_t)
      arr = filter(x -> x != e, elem_t)
      diff = norm(node[arr[1]] - node[arr[2]]) .^ 2
      edgelength[t, j] = diff
    end
  end
  max_indices = findmax(edgelength, dims = 2)[2]
  for t = 1:NT
    shifted = _shift_to_first(elem[t][:], T(max_indices[t][2]))
    elem[t] = shifted
  end
end

function _setup_markers_and_nodes!(
  node::AbstractVector{<:VectorValue},
  elem::AbstractVector{<:AbstractVector{T}},
  d2p::SparseMatrixCSC{T,T},
  dualedge::SparseMatrixCSC{T},
  NE::T,
  η_arr::AbstractArray,
  θ::AbstractFloat,
) where {T<:Integer}
  total_η = sum(η_arr)
  partial_η = 0
  sorted_η_idxs = sortperm(-η_arr)
  markers = zeros(T, NE)
  # Loop over global triangle indices
  for t in sorted_η_idxs
    if (partial_η >= θ * total_η)
      break
    end
    need_to_mark = true
    # Get triangle index with next largest error
    while (need_to_mark)
      # Base point
      base = d2p[elem[t][2], elem[t][3]]
      # Already marked
      if markers[base] > 0
        need_to_mark = false
      else
        # Get the estimator contribution for the current triangle
        partial_η = partial_η + η_arr[t]
        # Increase the number of nodes to add new midpoint
        N = size(node, 1) + 1
        # The markers of the current elements is this node
        markers[d2p[elem[t][2], elem[t][3]]] = N
        # Coordinates of new node
        elem2 = elem[t][2]
        elem3 = elem[t][3]
        midpoint = _get_midpoint(node[[elem2, elem3], :])
        node = push!(node, midpoint)
        t = dualedge[elem[t][3], elem[t][2]]
        # There is no dual edge here, go to next triangle index
        if t == 0
          need_to_mark = false
        end
      end
    end
  end
  node, markers
end

function _divide!(elem::AbstractVector, t::T, p::AbstractVector{T}) where {T<:Integer}
  new_row = [p[4], p[3], p[1]]
  update_row = [p[4], p[1], p[2]]
  push!(elem, new_row)
  elem[t] = update_row
  elem
end

function _bisect(
  d2p::SparseMatrixCSC{T,T},
  elem::AbstractVector,
  markers::AbstractVector{T},
  NT::T,
) where {T<:Integer}
  forest = Vector{BinaryNode}()
  for t in UnitRange{T}(1:NT)
    base = d2p[elem[t][2], elem[t][3]]
    if (markers[base] > 0)
      newnode = BinaryNode(t)
      p = vcat(elem[t][:], markers[base])
      elem = _divide!(elem, t, p)
      cur_size::T = size(elem, 1)
      # l and r are BinaryNodes for building forest
      l = leftchild(t, newnode)
      r = rightchild(cur_size, newnode)
      # left and right are point indices
      left = d2p[p[1], p[2]]
      right = d2p[p[3], p[1]]
      if (markers[right] > 0)
        leftchild(cur_size, r)
        # Need to handle edge case where both markers are positive.
        # In this case, we have to offset by 2 because the other branch
        # already took the next index.
        if markers[left] > 0
          rightchild(cur_size + 2, r)
        else
          rightchild(cur_size + 1, r)
        end
        elem = _divide!(elem, cur_size, [p[4], p[3], p[1], markers[right]])
      end
      if (markers[left] > 0)
        leftchild(t, l)
        rightchild(cur_size + 1, l)
        elem = _divide!(elem, t, [p[4], p[1], p[2], markers[left]])
      end
      push!(forest, newnode)
    else
      push!(forest, BinaryNode(t))
    end
  end
  elem, forest
end

function _build_edges(elem::AbstractVector{<:AbstractVector{T}}) where {T<:Integer}
  edge = zeros(T, 3 * length(elem), 2)
  for t = 1:length(elem)
    off = 3 * (t - 1)
    edge[off+1, :] = [elem[t][1] elem[t][2]]
    edge[off+2, :] = [elem[t][1] elem[t][3]]
    edge[off+3, :] = [elem[t][2] elem[t][3]]
  end
  edge = unique(sort!(edge, dims = 2), dims = 1)
end

function _build_directed_dualedge(
  elem::AbstractVector{<:AbstractVector{T}},
  N::T,
  NT::T,
) where {T<:Integer}
  dualedge = spzeros(T, N, N)
  for t = 1:NT
    dualedge[elem[t][1], elem[t][2]] = t
    dualedge[elem[t][2], elem[t][3]] = t
    dualedge[elem[t][3], elem[t][1]] = t
  end
  dualedge
end

function _dual_to_primal(edge::Matrix{T}, NE::T, N::T) where {T<:Integer}
  d2p = spzeros(T, T, N, N)
  for k = 1:NE
    i = edge[k, 1]
    j = edge[k, 2]
    d2p[i, j] = k
    d2p[j, i] = k
  end
  d2p
end

function _is_against_top(face::Table{<:Integer}, top::GridTopology, d::Integer)
  face_vec = sort.([face[i][:] for i = 1:size(face, 1)])
  face_top = sort.(get_faces(top, d, 0))
  issetequal_bitvec = issetequal(face_vec, face_top)
  all(issetequal_bitvec)
end

function _is_against_top(face::Matrix{<:Integer}, top::GridTopology, d::Integer)
  face_vec = sort.([face[i, :] for i = 1:size(face, 1)])
  face_top = sort.(get_faces(top, d, 0))
  issetequal_bitvec = issetequal(face_vec, face_top)
  all(issetequal_bitvec)
end

# For initial labeling
function _sort_ccw(cell_coords::AbstractVector{<:VectorValue})
  midpoint = _get_midpoint(cell_coords)
  offset_coords = cell_coords .- midpoint
  sortperm(offset_coords, by = v -> atan(v[2], v[1]))
end

# For initial labeling
function _sort_cell_node_ids_ccw!(
  cell_node_ids::AbstractVector{<:AbstractVector{T}},
  node_coords::AbstractVector{<:VectorValue},
) where {T<:Integer}
  for (i, cell) in enumerate(cell_node_ids)
    cell_coords = node_coords[cell]
    perm = _sort_ccw(cell_coords)
    cell_node_ids[i] = cell[perm]
  end
end

function _markers_to_node_to_bis_edge(markers, edge)
  node_to_bis_edge = Dict()
  for (i, marker) in enumerate(markers)
    # Edge was bisected
    if marker != 0
      node_to_bis_edge[marker] = edge[i, :]
    end
  end
  node_to_bis_edge
end

"""
Lowest level interface to the newest vertex bisection algorithm. This step
takes places after unpacking the Grid. It also uses the `should_sort` to
determine if an initial sorting by the longest edge (for nonuniform meshes) is
necessary.

`node_coords, cell_node_ids -> node_coords, cell_node_ids`

# Arguments

 -`node_coords::AbstractVector{<:VectorValue}`: The vector of d-dimensional
 nodal coordinates stored as `VectorValue`s

 -`cell_node_ids::AbstractVector{<:AbstractVector{T}}`: Contains the elements
 as ordered tuples of nodal indices.

 -`η_arr::AbstractArray`: The values of the estimator on each cell of the domain

 -`θ::AbstractArray`: Dörfler marking parameter: 0 = no refinement,
   1=uniform refinement.

"""
function newest_vertex_bisection(
  node_coords::AbstractVector{<:VectorValue},
  cell_node_ids::AbstractVector{<:AbstractVector{T}},
  η_arr::AbstractArray,
  θ::AbstractFloat,
) where {T<:Integer}
  # Number of nodes
  N::T = size(node_coords, 1)
  # Long Chen terminology
  elem = cell_node_ids
  # Number of cells (triangles in 2D)
  NT::T = size(elem, 1)
  @assert length(η_arr) == NT
  # Long Chen terminology
  edge = _build_edges(elem)
  NE::T = size(edge, 1)
  # Long Chen terminology
  dualedge = _build_directed_dualedge(elem, N, NT)
  d2p = _dual_to_primal(edge, NE, N)
  node_coords, markers =
    _setup_markers_and_nodes!(node_coords, elem, d2p, dualedge, NE, η_arr, θ)
  # For performance since we need to push! to elem
  elem = Vector{Vector}(elem)
  cell_node_ids, forest = _bisect(d2p, elem, markers, NT)
  node_to_bis_edge = _markers_to_node_to_bis_edge(markers, edge)
  node_coords, cell_node_ids, forest, node_to_bis_edge
end

"""
Middle level interface to the newest vertex bisection algorithm. This step
takes places after unpacking the DiscreteModel and performs sorting if
necessary. It maps

`Grid -> Grid, buffer`

# Arguments

 -`grid::Grid`: The current Grid.

 -`η_arr::AbstractArray`: The values of the estimator on each cell of the domain

 -`θ::AbstractArray`: Dörfler marking parameter: 0 = no refinement,
   1=uniform refinement.

"""
function newest_vertex_bisection(grid::Grid, η_arr::AbstractArray, θ::AbstractFloat)
  node_coords = get_node_coordinates(grid)
  # Need "un lazy" version for resize!
  node_coords = [v for v in node_coords]
  #@show cell_node_ids = get_cell_node_ids(grid)
  cell_node_ids = get_cell_node_ids(grid)
  # Convert from Table to Vector{Vector}
  cell_node_ids = [v for v in cell_node_ids]
  # Should always sort on the first iteration
  _sort_cell_node_ids_ccw!(cell_node_ids, node_coords)
  _sort_longest_edge!(cell_node_ids, node_coords)
  node_coords_ref, cell_node_ids_unsort, forest, node_to_bis_edge =
    newest_vertex_bisection(node_coords, cell_node_ids, η_arr, θ)
  reffes = get_reffes(grid)
  cell_types = get_cell_type(grid)
  cell_types = [c for c in cell_types]
  # TODO : Is this ok since we only handle triangular mesh for now?
  new_cell_types = fill(1, length(cell_node_ids_unsort) - length(cell_node_ids))
  append!(cell_types, new_cell_types)
  cell_node_ids_ref = Table([c for c in cell_node_ids_unsort])
  grid_ref_unsort = UnstructuredGrid(node_coords_ref, cell_node_ids_ref, reffes, cell_types)
  # Other buffer information can be added later
  buffer = (; grid_ref_unsort)
  sort!.(cell_node_ids_unsort)
  cell_node_ids_ref = Table([c for c in cell_node_ids_unsort])
  grid_ref = UnstructuredGrid(node_coords_ref, cell_node_ids_ref, reffes, cell_types)
  grid_ref, buffer, forest, node_to_bis_edge
end

"""
The newest vertex bisection algorithm provides a method of local refinement
without creating hanging nodes. For now, only 2D simplicial meshes are
supported.

This is the highest level version of the function, it maps

`DiscreteModel -> DiscreteModel, buffer`


# Arguments

 -`model::DiscreteModel`: The current DiscreteModel to be refined.

 -`η_arr::AbstractArray`: The values of the estimator on each cell of the domain
 i.e. one should have `length(η_arr) == num_cells(model)`

 -`θ::AbstractArray=1.0`: Dörfler marking parameter: 0 = no refinement,
   1=uniform refinement.

"""
function newest_vertex_bisection(
  model::DiscreteModel,
  η_arr::AbstractArray;
  θ = 1.0, # corresponds to uniform refinement
)
  @assert length(η_arr) == num_cells(model)
  grid = get_grid(model)
  topo = GridTopology(grid)
  grid_ref, buffer, forest, node_to_bis_edge = newest_vertex_bisection(grid, η_arr, θ)
  topo_ref = GridTopology(grid_ref)
  d_to_dface_to_olddim, d_to_dface_to_oldid = _create_d_to_dface_to_old(
    forest,
    topo,
    topo_ref,
    get_node_coordinates(grid_ref),
    node_to_bis_edge,
  )
  labels_ref = _propagate_labeling(model, d_to_dface_to_olddim, d_to_dface_to_oldid)
  DiscreteModel(grid_ref, topo_ref, labels_ref), buffer
end

"""
The newest vertex bisection algorithm provides a method of local refinement
without creating hanging nodes. For now, only 2D simplicial meshes are
supported.

This is the highest level version of the function, it maps

`DiscreteModel, buffer -> DiscreteModel, buffer`


# Arguments

 -`model::DiscreteModel`: The current DiscreteModel to be refined.

 -`buffer::NamedTuple`: The buffer providing all the itermediate information
 between refinements. For now just includes a Grid with orderings necessary
 for the algorithm.

 -`η_arr::AbstractArray`: The values of the estimator on each cell of the domain
 i.e. one should have `length(η_arr) == num_cells(model)`

 -`θ::AbstractArray=1.0`: Dörfler marking parameter: 0 = no refinement,
   1=uniform refinement.

"""
function newest_vertex_bisection(
  model::DiscreteModel,
  buffer::NamedTuple,
  η_arr::AbstractArray;
  θ = 1.0, # corresponds to uniform refinement
)
  @assert length(η_arr) == num_cells(model)
  grid = get_grid(model)
  grid_ref_unsort = buffer.grid_ref_unsort
  topo = GridTopology(grid)
  grid_ref, buffer, forest, node_to_bis_edge =
    newest_vertex_bisection(grid_ref_unsort, η_arr, θ)
  topo_ref = GridTopology(grid_ref)
  d_to_dface_to_olddim, d_to_dface_to_oldid = _create_d_to_dface_to_old(
    forest,
    topo,
    topo_ref,
    get_node_coordinates(grid_ref),
    node_to_bis_edge,
  )
  labels_ref = _propagate_labeling(model, d_to_dface_to_olddim, d_to_dface_to_oldid)
  model_ref = DiscreteModel(grid_ref, topo_ref, labels_ref)
  model_ref, buffer
end
