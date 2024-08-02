"""
Note on RefinementRules and Orientation of the refined grids:

  In order to guarantee that the refined grid is Oriented (which is something
  we want for div- and curl-conforming discretisations), we need to guarantee
  for simplices that each fine cell has it's vertex gids sorted in increasing
  order.

  In the case of refined meshes, this has an additional constraint: the sorting
  of the gids CANNOT be done after the refinement. This would make the glue
  inconsistent. This is an attempt to explain why:
  If we change the order of the gids of a fine cell, we are basically applying a
  rotation(+symmetry) to the reference space of that cell. The mesh itself will be
  aware of this rotation, but the RefinementRule will not. So the reference space
  corresponding to that fine cell within the RefinementRule will NOT be rotated
  accordingly. This creates an inconsistency, and the fine-to-coarse/coarse-to-fine
  mesh transfers will therefore be broken (the glue is no longer valid).

  How to fix this in the future? The glue should probably keep a record of the
  cell-wise rotations. This could be an auxiliary cell_map that maps the reference
  space of the fine cell within the fine mesh to the reference space of the fine cell
  within the refinement rule.
  For instance:
    - Φ: Cell_map given by the RefinementRule, going from the fine reference space
         to the coarse reference space.
    - β: Cell_map encoding the rotation from the reference space of the fine mesh to
         the reference space of the RefinementRule.
  then we would have
      X_coarse = Φ∘β(X_fine)
"""

abstract type EdgeBasedCoarsening end

struct NVBCoarsening{T} <: EdgeBasedCoarsening
    cell_to_longest_edge_lid::Vector{T}
    cell_to_longest_edge_gid::Vector{T}
end

function NVBCoarsening(model::DiscreteModel{Dc,Dp}) where {Dc,Dp}
    topo = model.grid_topology
    @check all(p -> p == TRI, get_polytopes(topo))
    c2e_map = get_faces(topo, Dc, 1)
    e2n_map = get_faces(topo, 1, 0)
    node_coords = get_node_coordinates(model)
    longest_edge_lids, longest_edge_gids =
        _get_longest_edge_ids(c2e_map, e2n_map, node_coords)
    NVBCoarsening(longest_edge_lids, longest_edge_gids)
end

NVBCoarsening(model::AdaptedDiscreteModel) = NVBRefinement(model.model)

function coarsen(
    method::EdgeBasedCoarsening,
    model::UnstructuredDiscreteModel{Dc,Dp},
    glue::AdaptivityGlue;
    cells_to_coarsen = nothing,
) where {Dc,Dp}
    # Create new mode
    rrules, faces_list =
        setup_edge_based_crules(method, model.grid_topology, glue, cells_to_coarsen)
    topo = _coarsen_unstructured_topology(method, model.grid_topology, glue)
    reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
    grid   = UnstructuredGrid(get_vertex_coordinates(topo),get_faces(topo,Dc,0),reffes,get_cell_type(topo),OrientationStyle(topo))
    labels = FaceLabeling(topo)
    ref_model = UnstructuredDiscreteModel(grid,topo,labels)
    #ref_model = model
    #topo = get_grid_topology(model)
    ## Create ref glue
    glue = _get_refinement_glue(topo, ref_model.grid_topology, rrules)
    return AdaptedDiscreteModel(ref_model, model, glue)
end

function _get_refinement_glue(ftopo ::UnstructuredGridTopology{Dc},
                              ctopo ::UnstructuredGridTopology{Dc},
                              rrules::AbstractVector{<:RefinementRule}) where {Dc}
  nC_old = num_faces(ctopo,Dc)
  nC_new = num_faces(ftopo,Dc)

  f2c_cell_map      = Vector{Int}(undef,nC_new)
  fcell_to_child_id = Vector{Int}(undef,nC_new)

  k = 1
  for iC = 1:nC_old
    rr = rrules[iC]
    range = k:k+num_subcells(rr)-1
    f2c_cell_map[range] .= iC
    fcell_to_child_id[range] .= collect(1:num_subcells(rr))
    k += num_subcells(rr)
  end

  f2c_faces_map = [(d==Dc) ? f2c_cell_map : Int[] for d in 0:Dc]
  return AdaptivityGlue(f2c_faces_map,fcell_to_child_id,rrules)
end

function _is_newest_vertex(n, e2n_map_cache, e2n_map, longest_edge_gid)
    longest_edge_nodes = getindex!(e2n_map_cache, e2n_map, longest_edge_gid)
    #@show longest_coords
    return n ∉ longest_edge_nodes
end

function _identify_good_nodes_to_coarsen(
    n2c_map_cache,
    n2c_map,
    topo,
    nN,
    c_to_longest_edge_gid,
)
    e2n_map = get_faces(topo, 1, 0)
    e2n_map_cache = array_cache(e2n_map)
    good_nodes = typeof(nN)[]
    for n = 1:nN
        patch_cells = getindex!(n2c_map_cache, n2c_map, n)
        if length(patch_cells) == 2 || length(patch_cells) == 4
            is_newest = Bool[]
            for c in patch_cells
                #@show c
                #@show coords[getindex!(c2n_map_cache, c2n_map, c)]
                push!(
                    is_newest,
                    _is_newest_vertex(n, e2n_map_cache, e2n_map, c_to_longest_edge_gid[c]),
                )
            end
            if all(is_newest)
                push!(good_nodes, n)
            end
        end
    end
    return good_nodes
end

function _combine_brothers(
    brothers::Vector{<:Pair},
    e2n_map,
    e2n_map_cache,
    c_to_longest_edge_gid,
    coords,
)
    new_coarse_nodes = Vector{Vector{eltype(c_to_longest_edge_gid)}}()
    for (c1, c2) in brothers
        nodes = Set{eltype(c_to_longest_edge_gid)}()
        #@show (c1, c2)
        e1 = c_to_longest_edge_gid[c1]
        e2 = c_to_longest_edge_gid[c2]
        nodes_c1 = getindex!(e2n_map_cache, e2n_map, e1)
        union!(nodes, nodes_c1)
        #@show nodes_c1
        nodes_c2 = getindex!(e2n_map_cache, e2n_map, e2)
        union!(nodes, nodes_c2)
        #nodes_c1 = getindex!(c2n_map_cache, c2n_map,  c1)
        #nodes_c2 = getindex!(c2n_map_cache, c2n_map,  c2)
        #@show nodes_c2
        coarse_nodes = collect(nodes)
        push!(new_coarse_nodes, coarse_nodes)
        my_coords = coords[coarse_nodes]
    end
    return new_coarse_nodes
end

function _get_brother_pairs(good_nodes, n2c_map_cache, n2c_map, c2e_map_cache, c2e_map, cell_to_parent_gid)
    all_brother_pairs = Dict{
        eltype(good_nodes),
        Vector{Pair{eltype(cell_to_parent_gid),eltype(cell_to_parent_gid)}},
    }()
    intersect_empty(a, b) = isempty(intersect(a, b))
    for n in good_nodes
        patch_cells = getindex!(n2c_map_cache, n2c_map, n)
        #@show n
        parent_cell_gids = eltype(good_nodes)[]
        for c in patch_cells
            #@show cell_to_parent_gid[c]
            push!(parent_cell_gids, cell_to_parent_gid[c])
        end
        #@show unique!(parent_cell_gids)
        if length(patch_cells) == 2
            brother_pairs = [Pair(patch_cells[1], patch_cells[2])]
        elseif length(patch_cells) == 4
            edges1 = getindex!(c2e_map_cache, c2e_map, patch_cells[1]) |> copy
            edges2 = getindex!(c2e_map_cache, c2e_map, patch_cells[2]) |> copy
            edges3 = getindex!(c2e_map_cache, c2e_map, patch_cells[3]) |> copy
            if !intersect_empty(edges1, edges2)#!isempty(intersect(edges1, edges2))
                brother_pairs =
                [Pair(patch_cells[1], patch_cells[2]), Pair(patch_cells[3], patch_cells[4])]
            elseif !intersect_empty(edges1, edges3)
                brother_pairs =
                [Pair(patch_cells[1], patch_cells[3]), Pair(patch_cells[2], patch_cells[4])]
            else
                brother_pairs =
                [Pair(patch_cells[1], patch_cells[4]), Pair(patch_cells[2], patch_cells[3])]
            end
        else
            error("Number of edges not equal to 2 or 4 :)")
        end
        all_brother_pairs[n] = brother_pairs
    end
    return all_brother_pairs
end

function setup_edge_based_crules(
    method::NVBCoarsening,
    topo::UnstructuredGridTopology{Dc},
    glue::AdaptivityGlue,
    cells_to_coarsen::AbstractArray{<:Integer},
) where {Dc}
    nC = num_faces(topo, Dc)
    return [WhiteRefinementRule(TRI) for _ = 1:nC], nothing
end

function _coarsen_unstructured_topology(
    method::EdgeBasedCoarsening,
    topo::UnstructuredGridTopology{Dc},
    glue::AdaptivityGlue,
) where {Dc}
    nC = num_faces(topo, Dc)
    nE = num_faces(topo, 1)
    nN = num_faces(topo, 0)
    n2c_map = get_faces(topo, 0, Dc)
    n2c_map_cache = array_cache(n2c_map)
    c2n_map = get_faces(topo, Dc, 0)
    c2n_map_cache = array_cache(n2c_map)
    c2e_map = get_faces(topo, Dc, 1)
    e2n_map = get_faces(topo, 1, 0)
    e2n_map_cache = array_cache(e2n_map)
    c2e_map_cache = array_cache(c2e_map)
    e2c_map = get_faces(topo, 1, Dc)
    polys = topo.polytopes
    cell_types = topo.cell_type
    cell_color = copy(cell_types) # WHITE
    coords = get_vertex_coordinates(topo)
    c_to_longest_edge_gid = method.cell_to_longest_edge_gid
    good_nodes = _identify_good_nodes_to_coarsen(
        n2c_map_cache,
        n2c_map,
        topo,
        nN,
        c_to_longest_edge_gid,
    )
    coords_new = setdiff(coords, coords[good_nodes]) 
    cell_to_parent_gid = glue.n2o_faces_map[3]
    all_brother_pairs =
        _get_brother_pairs(good_nodes, n2c_map_cache, n2c_map, c2e_map_cache, c2e_map, cell_to_parent_gid)
    cells_to_remove = eltype(cell_to_parent_gid)[]
    for (good_node, pairs) in all_brother_pairs
        for pair in pairs
            push!(cells_to_remove, pair.first)
            push!(cells_to_remove, pair.second)
        end
    end
    new_cells = []
    for brothers in values(all_brother_pairs)
        union!(new_cells, _combine_brothers(brothers, e2n_map, e2n_map_cache, c_to_longest_edge_gid, coords))
    end
    # Setup c2n_map_new
    num_new_cells = length(c2n_map) - length(cells_to_remove)
    c2n_map_new = Vector{Vector{eltype(cell_to_parent_gid)}}(undef,num_new_cells)
    c_new = 0
    for c = 1:nC
        # If this cell should not be removed
        if c ∉ cells_to_remove
            nodes = copy(getindex!(c2n_map_cache, c2n_map, c))
            c_new += 1
           c2n_map_new[c_new] = nodes
       end
    end
    union!(c2n_map_new, new_cells)
    for nodes in c2n_map_new
        for i in eachindex(nodes)
            offset = filter(n -> nodes[i] > n, good_nodes) |> length
            nodes[i] -= offset
        end
    end
    c2n_map_new = Table(c2n_map_new)
    #coords_new = get_new_coordinates_from_faces(topo, faces_list)
    #c2n_map_new = get_refined_cell_to_vertex_map(topo, rrules, faces_list)
    #polys_new, cell_type_new = _get_cell_polytopes(rrules)

    ## We can guarantee the new topology is oriented if
    ##   1 - the old topology was oriented
    ##   2 - we have a single type of polytope (i.e new topo is not mixed)
    ##orientation = NonOriented()
    #orientation = NonOriented()
    new_cell_type = ones(Int8, length(c2n_map_new))
    return UnstructuredGridTopology(
        coords_new,
        c2n_map_new,
        new_cell_type,
        topo.polytopes,
        Oriented(),
    )
end
