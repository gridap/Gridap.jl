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

abstract type EdgeBasedRefinement <: AdaptivityMethod end


struct NVBRefinement{T} <: EdgeBasedRefinement
  cell_to_longest_edge_lid::Vector{T}
  cell_to_longest_edge_gid::Vector{T}
end

_calculate_edge_length(node_coords, nodes) = norm(node_coords[nodes[1]] - node_coords[nodes[2]])

# Create Vectors containing the longest edge gids and lids
function _get_longest_edge_ids(c2e_map, e2n_map, node_coords)
  nC = length(c2e_map)
  e2n_map_cache = array_cache(e2n_map)
  c2e_map_cache = array_cache(c2e_map)
  longest_edge_gids = Vector{eltype(eltype(c2e_map))}(undef, nC)
  longest_edge_lids = Vector{eltype(eltype(c2e_map))}(undef, nC)
  # Loop over cells
  for c in 1:length(c2e_map)
    es = getindex!(c2e_map_cache, c2e_map, c)
    e_length_max = 0.0
    e_gid_max = -1 # Longest edge index must be found
    e_lid_max = -1 # Longest edge index must be found
    for (e_lid, e_gid) in enumerate(es)
      ns = getindex!(e2n_map_cache, e2n_map, e_gid)
      e_length = _calculate_edge_length(node_coords, ns)
      if e_length > e_length_max
        e_length_max = e_length
        e_gid_max = e_gid
        e_lid_max = e_lid
      end
    end
    longest_edge_gids[c] = e_gid_max
    longest_edge_lids[c] = e_lid_max
  end
  return longest_edge_lids, longest_edge_gids
end

function NVBRefinement(model::DiscreteModel{Dc,Dp}) where {Dc, Dp}
  topo = model.grid_topology
  @check all(p->p==TRI,get_polytopes(topo))
  c2e_map     = get_faces(topo,Dc,1)
  e2n_map     = get_faces(topo,1 ,0)
  node_coords = get_node_coordinates(model)
  longest_edge_lids, longest_edge_gids = _get_longest_edge_ids(c2e_map, e2n_map, node_coords)
  NVBRefinement(longest_edge_lids, longest_edge_gids)
end

NVBRefinement(model::AdaptedDiscreteModel) = NVBRefinement(model.model)

struct RedGreenRefinement <: EdgeBasedRefinement end

function refine(method::EdgeBasedRefinement,model::UnstructuredDiscreteModel{Dc,Dp};cells_to_refine=nothing) where {Dc,Dp}
  # cells_to_refine can be
  #    a) nothing -> All cells get refined
  #    b) AbstractArray{<:Bool} of size num_cells(model)
  #            -> Only cells such that cells_to_refine[iC] == true get refined
  #    c) AbstractArray{<:Integer}
  #            -> Cells for which gid ∈ cells_to_refine get refined

  # Create new model
  rrules, faces_list = setup_edge_based_rrules(method, model.grid_topology,cells_to_refine)
  topo   = _refine_unstructured_topology(model.grid_topology,rrules,faces_list)
  reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  grid   = UnstructuredGrid(get_vertex_coordinates(topo),get_faces(topo,Dc,0),reffes,get_cell_type(topo),OrientationStyle(topo))
  labels = FaceLabeling(topo)
  ref_model = UnstructuredDiscreteModel(grid,topo,labels)
  ## Create ref glue
  glue = _get_refinement_glue(topo,model.grid_topology,rrules)

  return AdaptedDiscreteModel(ref_model,model,glue)
end

function _refine_unstructured_topology(topo::UnstructuredGridTopology{Dc},
                                      rrules::AbstractVector{<:RefinementRule},
                                      faces_list::Tuple) where {Dc}
  coords_new  = get_new_coordinates_from_faces(topo,faces_list)
  c2n_map_new = get_refined_cell_to_vertex_map(topo,rrules,faces_list)
  polys_new, cell_type_new = _get_cell_polytopes(rrules)

  # We can guarantee the new topology is oriented if
  #   1 - the old topology was oriented
  #   2 - we have a single type of polytope (i.e new topo is not mixed)
  #orientation = NonOriented()
  orientation = NonOriented()

  return UnstructuredGridTopology(coords_new,c2n_map_new,cell_type_new,polys_new,orientation)
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

"""
Given an UnstructuredTopology and a list of cells_to_refine, provides an edge-based coloring
for the topology:

  - If a cell is completely refined, it is colored RED
  - If a cell touches more than one refined (RED) cell, it becomes RED.
  - If a cell touches a single refined (RED) cell, it is colored GREEN.

The method returns a vector of Red/Green refinement rules, as well the list of
vertices, edges and cells which correspond to new vertices in the refined topology.
"""
function setup_edge_based_rrules(::EdgeBasedRefinement, topo::UnstructuredGridTopology{Dc},::Nothing) where Dc
  rrules     = lazy_map(RedRefinementRule,CompressedArray(topo.polytopes,topo.cell_type))
  faces_list = _redgreen_refined_faces_list(topo,rrules,[true])
  return rrules, faces_list
end

function setup_edge_based_rrules(::EdgeBasedRefinement,topo::UnstructuredGridTopology{Dc},cells_to_refine::AbstractArray{<:Bool}) where Dc
  return setup_edge_based_rrules(topo,findall(cells_to_refine))
end

_shift_to_first(v::AbstractVector{T}, i::T) where {T<:Integer} = circshift(v, -(i - 1))

function setup_edge_based_rrules(method::NVBRefinement, topo::UnstructuredGridTopology{Dc},cells_to_refine::AbstractArray{<:Integer}) where Dc
  nE = num_faces(topo,1)
  c2e_map       = get_faces(topo,Dc,1)
  c2e_map_cache = array_cache(c2e_map)
  e2c_map       = get_faces(topo,1,Dc)
  polys       = topo.polytopes
  cell_types  = topo.cell_type
  cell_color  = copy(cell_types) # WHITE
  # Hardcoded for TRI
  WHITE, GREEN, BLUE, BLUE_DOUBLE = Int8(1), Int8(2), Int8(5), Int(11)
  # Create the inverse mapping from long/short edge pairs to unique indices
  BLUE_dict = Dict((1,2) => 1, (1, 3) => 2, (2, 1) => 3, (2, 3) => 4, (3, 1) => 5, (3, 2) => 6)
  green_offset = 3
  blue_offset = 3*2
  blue_double_offset = 3
  c_to_longest_edge_gid = method.cell_to_longest_edge_gid
  c_to_longest_edge_lid = method.cell_to_longest_edge_lid
  is_refined = falses(nE)
  # Loop over cells and mark edges to refine i.e. is_refined
  # The reason to not loop directly on c is that we need to change c within
  # a given iteration of the for loop
  for i in 1:length(cells_to_refine)
    c = cells_to_refine[i]
    e_longest = c_to_longest_edge_gid[c]
    # Has to terminate because an edge is marked each iteration or we skip an
    # iteration due to a boundary cell
    while !is_refined[e_longest]
      is_refined[e_longest] = true
      c_nbor_lid = findfirst(c′ -> c′ != c, e2c_map[e_longest])
      if isnothing(c_nbor_lid) # We've reach the boundary
        continue
      else
        # Get the longest edge of the neighbor
        c_nbor_gid = e2c_map[e_longest][c_nbor_lid]
        e_longest = c_to_longest_edge_gid[c_nbor_gid]
        # Set the current cell gid to that of the neighbor
        c = c_nbor_gid
      end
    end
  end
  # Loop over cells and refine based on marked edges
  for c in 1:length(c2e_map)
    c_edges = getindex!(c2e_map_cache, c2e_map, c)
    refined_edge_lids = findall(is_refined[c_edges])
    # GREEN refinement because only one edge should be bisected
    if length(refined_edge_lids) == 1
      ref_edge = refined_edge_lids[1]
      cell_color[c] = GREEN + Int8(ref_edge-1)
      # BLUE refinement: two bisected
    elseif length(refined_edge_lids) == 2
      long_ref_edge_lid = c_to_longest_edge_lid[c]
      short_ref_edge_lid = setdiff(refined_edge_lids, long_ref_edge_lid)[1]
      blue_idx = BLUE_dict[(long_ref_edge_lid, short_ref_edge_lid)]
      cell_color[c] = BLUE + Int8(blue_idx - 1)
      # DOUBLE BLUE refinement: three bisected edges (somewhat rare)
    elseif length(refined_edge_lids) == 3
      long_ref_edge_lid = c_to_longest_edge_lid[c]
      cell_color[c] = BLUE_DOUBLE + Int(long_ref_edge_lid - 1)
    end
  end
  # Setup color_rrules for the CompressedArray
  num_rr =  1 + green_offset + blue_offset + blue_double_offset # GREEN+BLUE+DOUBLE_BLUE
  T = typeof(WhiteRefinementRule(first(polys))) # Needed to make eltype(color_rrules) concrete
  color_rrules = Vector{T}(undef,num_rr)
  p = polys[1] # Hardcoded for TRI
  color_rrules[WHITE] = WhiteRefinementRule(p)
  for e in 1:num_faces(p,1)
    color_rrules[GREEN+e-1] = GreenRefinementRule(p,e)
  end
  # blue_idx corresponds to offset into BLUE_dict
  blue_idx = 0
  for e1 in 1:num_faces(p, 1)
    for e2 in 1:num_faces(p, 1)
      if e1 != e2
        color_rrules[BLUE + blue_idx] = BlueRefinementRule(p, e1, e2)
        blue_idx += 1
      end
    end
  end
  for e in 1:num_faces(p, 1)
    color_rrules[BLUE_DOUBLE + e - 1] = BlueDoubleRefinementRule(p, e)
  end
  rrules = CompressedArray(color_rrules,cell_color)
  face_list = _redgreen_refined_faces_list(topo, rrules, is_refined)
  return rrules, face_list
end

function setup_edge_based_rrules(::RedGreenRefinement, topo::UnstructuredGridTopology{Dc},cells_to_refine::AbstractArray{<:Integer}) where Dc
  nC = num_cells(topo)
  nE = num_faces(topo,1)
  c2e_map     = get_faces(topo,Dc,1)
  e2c_map     = get_faces(topo,1,Dc)
  polys       = topo.polytopes
  cell_types  = topo.cell_type
  nP          = length(polys)

  # Color pointers/offsets
  WHITE, RED, GREEN = Int8(1), Int8(1+nP), Int8(1+nP+nP)
  green_offsets = counts_to_ptrs(map(p->num_faces(p,1),polys))

  cell_color  = copy(cell_types) # WHITE
  is_red      = fill(false,nC)
  is_refined  = fill(false,nE)

  # Flag initial red cells and edges
  for c in cells_to_refine
    cell_color[c] = RED + Int8(cell_types[c]-1)
    is_red[c] = true
    is_refined[c2e_map[c]] .= true
  end

  # Propagate red/green flags
  # Queue invariant: Cells in queue are RED, every RED cell is only touched once
  q = Queue{Int}()
  map(c->enqueue!(q,c),cells_to_refine)
  while !isempty(q)
    c = dequeue!(q)
    c_edges = c2e_map[c]

    # For each non-red neighboring cell
    nbors = map(first,filter(n->n!=[] && !is_red[n[1]],map(e->filter(x->x!=c,e2c_map[e]),c_edges)))
    for nbor in nbors
      nbor_edges = c2e_map[nbor]
      nbor_ref_edges = findall(is_refined[nbor_edges])
      if length(nbor_ref_edges) == 1
        p = cell_types[nbor]
        ref_edge = nbor_ref_edges[1]
        cell_color[nbor] = GREEN + Int8(green_offsets[p]-1+ref_edge-1)
      else # > 1, can't be 0 (queue invariant)
        cell_color[nbor] = RED + Int8(cell_types[nbor]-1)
        is_red[nbor] = true
        is_refined[nbor_edges] .= true
        enqueue!(q,nbor)
      end
    end
  end

  # Create RefinementRules
  num_rr = nP + nP + green_offsets[nP+1]-1 # WHITE+RED+GREEN
  T = typeof(WhiteRefinementRule(first(polys))) # Needed to make eltype(color_rrules) concrete
  color_rrules = Vector{T}(undef,num_rr)
  for (k,p) in enumerate(polys)
    color_rrules[WHITE+k-1] = WhiteRefinementRule(p)
    color_rrules[RED+k-1]   = RedRefinementRule(p)
    for e in 1:num_faces(p,1)
      color_rrules[GREEN+green_offsets[k]-1+e-1] = GreenRefinementRule(p,e)
    end
  end
  rrules = CompressedArray(color_rrules,cell_color)

  # Create list of refined faces
  faces_list = _redgreen_refined_faces_list(topo,rrules,is_refined)

  return rrules, faces_list
end

function _redgreen_refined_faces_list(topo::UnstructuredGridTopology{2},rrules,edge_is_refined)
  ref_nodes = 1:num_faces(topo,0)
  ref_edges = all(edge_is_refined) ? (1:num_faces(topo,1)) : findall(edge_is_refined)
  ref_cells = findall(lazy_map(_has_interior_point,rrules))
  return (ref_nodes,ref_edges,ref_cells)
end

function get_new_coordinates_from_faces(p::Union{Polytope{D},GridTopology{D}},faces_list::Tuple) where {D}
  @check length(faces_list) == D+1

  nN_new     = sum(x->length(x),faces_list)
  coords_old = get_vertex_coordinates(p)
  coords_new = Vector{eltype(coords_old)}(undef,nN_new)

  n = 1
  for (d,dfaces) in enumerate(faces_list)
    if length(dfaces) > 0
      nf = length(dfaces)
      d2n_map = get_faces(p,d-1,0)
      coords_new[n:n+nf-1] .= map(f -> sum(coords_old[d2n_map[f]])/length(d2n_map[f]), dfaces)
      n += nf
    end
  end

  return coords_new
end

function get_refined_cell_to_vertex_map(topo::UnstructuredGridTopology{2},
                                        rrules::AbstractVector{<:RefinementRule},
                                        faces_list::Tuple)
  @notimplementedif !all(map(p -> p ∈ [QUAD,TRI],topo.polytopes))
  nN_old,nE_old,nC_old  = num_faces(topo,0),num_faces(topo,1),num_faces(topo,2)
  c2n_map = get_faces(topo,2,0)
  c2e_map = get_faces(topo,2,1)
  ref_edges = faces_list[2]

  # Allocate map ptr and data arrays
  nC_new    = sum(rr -> num_subcells(rr), rrules)
  nData_new = sum(rr -> sum(ids->length(ids),rr.ref_grid.grid.cell_node_ids), rrules)

  ptrs_new  = Vector{Int}(undef,nC_new+1)
  data_new  = Vector{Int}(undef,nData_new)

  cell_is_refined = lazy_map(rr->isa(RefinementRuleType(rr),RedRefinement),rrules)
  if all(cell_is_refined)
    edges_reindexing = 1:nE_old
  else
    edges_reindexing = lazy_map(Reindex(find_inverse_index_map(ref_edges,nE_old)),1:nE_old)
  end

  k = 1
  ptrs_new[1] = 1
  C = nN_old + length(ref_edges) + 1
  for iC = 1:nC_old
    rr = rrules[iC]

    # New Node gids from old N,E,C lids
    N = c2n_map[iC]
    E = edges_reindexing[c2e_map[iC]] .+ nN_old
    sub_conn = get_relabeled_connectivity(rr,(N,E,[C]))

    nChild = length(sub_conn)
    ptrs_new[k:k+nChild] .=  sub_conn.ptrs .+ (ptrs_new[k] - 1)
    data_new[ptrs_new[k]:ptrs_new[k+nChild]-1] .= sub_conn.data
    k = k+nChild

    _has_interior_point(rr) && (C += 1)
  end

  return Table(data_new,ptrs_new)
end

###############################################################################
# EdgeBasedRefinementRules

abstract type EdgeBasedRefinementRule <: RefinementRuleType end
struct RedRefinement <: EdgeBasedRefinementRule end
struct BlueRefinement{P, Q} <: EdgeBasedRefinementRule end
struct BlueDoubleRefinement{P} <: EdgeBasedRefinementRule end
struct GreenRefinement{P} <: EdgeBasedRefinementRule end

_has_interior_point(rr::RefinementRule) = _has_interior_point(rr,RefinementRuleType(rr))
_has_interior_point(rr::RefinementRule,::RefinementRuleType) = false

"""
RefinementRule representing a non-refined cell.
"""
function WhiteRefinementRule(p::Polytope)
  ref_grid = UnstructuredDiscreteModel(UnstructuredGrid(LagrangianRefFE(Float64,p,1)))
  return RefinementRule(WithoutRefinement(),p,ref_grid)
end

"""
Edge-based RefinementRule where a new vertex is added to
each edge of the original Polytope.
"""
function RedRefinementRule(p::Polytope{2})
  @notimplementedif (p ∉ [TRI,QUAD])

  if p == TRI
    faces_list = ([1,2,3],[1,2,3],[])
  elseif p == QUAD
    faces_list = ([1,2,3,4],[1,2,3,4],[1])
  end
  coords = get_new_coordinates_from_faces(p,faces_list)

  polys, cell_types, conn = _get_red_refined_connectivity(p)
  reffes = map(x->LagrangianRefFE(Float64,x,1),polys)

  ref_grid = UnstructuredDiscreteModel(UnstructuredGrid(coords,conn,reffes,cell_types))
  return RefinementRule(RedRefinement(),p,ref_grid)
end

function _get_red_refined_connectivity(p::Polytope{2})
  if p == TRI
    polys     = [TRI]
    cell_type = [1,1,1,1]
    conn_data = [1,4,5,
                2,4,6,
                3,5,6,
                4,5,6]
    conn_ptrs = [1,4,7,10,13]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  elseif p == QUAD
    polys     = [QUAD]
    cell_type = [1,1,1,1]
    conn_data = [1,5,7,9,
                5,2,9,8,
                7,9,3,6,
                9,8,6,4]
    conn_ptrs = [1,5,9,13,17]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  end
  @notimplemented
end

function _has_interior_point(rr::RefinementRule,::RedRefinement)
  p = get_polytope(rr)
  if p == QUAD
    return true
  end
  return false
end

"""
Edge-based RefinementRule where a new vertex is added to
a single edge of the original Polytope.
"""
function GreenRefinementRule(p::Polytope{2},ref_edge::Integer)
  @notimplementedif (p ∉ [TRI,QUAD])

  faces_list = (collect(1:num_faces(p,0)),[ref_edge],[])
  coords = get_new_coordinates_from_faces(p,faces_list)

  polys, cell_types, conn = _get_green_refined_connectivity(p,ref_edge)
  reffes = map(x->LagrangianRefFE(Float64,x,1),polys)

  ref_grid = UnstructuredDiscreteModel(UnstructuredGrid(coords,conn,reffes,cell_types))
  return RefinementRule(GreenRefinement{ref_edge}(),p,ref_grid)
end

function _get_green_refined_connectivity(p::Polytope{2},ref_edge)
  # Note: Sorting is necessary in order to guarantee that the gids
  #       of the refined mesh are sorted (and therefore that the fine
  #       grid is Oriented). See the note at top of the file.
  P = _get_green_vertex_permutation(p,ref_edge)
  if p == TRI
    polys     = [TRI]
    cell_type = [1,1]
    conn_data = [sort([4,P[1],P[2]])...,
                 sort([4,P[3],P[1]])...]
    conn_ptrs = [1,4,7]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  elseif p == QUAD
    polys     = [TRI]
    cell_type = [1,1,1]
    conn_data = [sort([P[1],P[2],5])...,
                 sort([P[3],P[1],5])...,
                 sort([P[4],5,P[2]])...]
    conn_ptrs = [1,4,7,10]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  end
  @notimplemented
end

function _get_green_vertex_permutation(p::Polytope{2},ref_edge::Integer)
  if p == TRI
    perm = circshift([1,2,3],ref_edge-3)
  elseif p == QUAD
    # Vertices and edges are not circularly labeled in QUAD, so
    # we have to be inventive:
    nodes = [1,2,4,3]  # Vertices in circular order
    nperm = [2,0,-1,1] # Number of perms depending on ref edge
    perm  = circshift(nodes,nperm[ref_edge])[nodes]
  else
    @notimplemented
  end
  return perm
end

"""
Edge-based RefinementRule where two vertex edges are added,
where the long edge is bisected first
"""
function BlueRefinementRule(p::Polytope{2}, long_ref_edge::Integer, short_ref_edge::Integer)
  @notimplementedif (p ∉ [TRI])

  faces_list = (collect(1:num_faces(p,0)),[long_ref_edge, short_ref_edge],[])
  coords = get_new_coordinates_from_faces(p,faces_list)
  polys, cell_types, conn = _get_blue_refined_connectivity(p,long_ref_edge, short_ref_edge)
  reffes = map(x->LagrangianRefFE(Float64,x,1),polys)

  ref_grid = UnstructuredDiscreteModel(UnstructuredGrid(coords,conn,reffes,cell_types))
  return RefinementRule(BlueRefinement{long_ref_edge, short_ref_edge}(),p,ref_grid)
end

function _get_blue_refined_connectivity(p::Polytope{2}, long_ref_edge, short_ref_edge)
  # Only implemented for triangles
  if p == TRI
    polys     = [TRI]
    cell_type = [1, 1, 1]
    unmarked_edge = setdiff([1,2,3], [long_ref_edge, short_ref_edge])[1]
    # Correspondence between an edge and the vertex opposite it
    edge_to_opposite_vertex = Dict(1 => 3, 2 => 2, 3 => 1)
    e2on = edge_to_opposite_vertex
    conn_data = [e2on[long_ref_edge], 4, 5,
                 e2on[unmarked_edge], 4, 5,
                 sort([e2on[long_ref_edge], e2on[short_ref_edge], 4])...]
    conn_ptrs = [1, 4, 7, 10]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  end
  @notimplemented
end

function _get_blue_vertex_permutation(p::Polytope{2},ref_edge::Integer)
  if p == TRI
    perm = circshift([1,2,3],ref_edge-3)
  else
    @notimplemented
  end
  return perm
end

"""
Edge-based RefinementRule where three vertex edges are added,
where the long edge is bisected first
"""
function BlueDoubleRefinementRule(p::Polytope{2}, long_ref_edge::Integer)
  @notimplementedif (p ∉ [TRI])
  otherfaces = setdiff(collect(1:num_faces(p,0)), long_ref_edge) |> sort
  faces_list = (collect(1:num_faces(p,0)),[otherfaces..., long_ref_edge],[])
  coords = get_new_coordinates_from_faces(p,faces_list)
  polys, cell_types, conn = _get_blue_double_refined_connectivity(p,long_ref_edge)
  reffes = map(x->LagrangianRefFE(Float64,x,1),polys)
  ref_grid = UnstructuredDiscreteModel(UnstructuredGrid(coords,conn,reffes,cell_types))
  return RefinementRule(BlueDoubleRefinement{long_ref_edge}(),p,ref_grid)
end

function _get_blue_double_refined_connectivity(p::Polytope{2}, long_ref_edge)
  if p == TRI
    polys     = [TRI]
    cell_type = [1, 1, 1, 1]
    # Hardcoded since there are only three cases
    if long_ref_edge == 3
      conn_data = [3, 5, 6,
                   1, 5, 6,
                   1, 4, 6,
                   2, 4, 6]
    elseif long_ref_edge == 2
      conn_data = [3, 5, 6,
                   2, 5, 6,
                   2, 4, 6,
                   1, 4, 6]
    elseif long_ref_edge == 1
      conn_data = [1, 4, 6,
                   3, 4, 6,
                   3, 5, 6,
                   2, 5, 6]
    end
    conn_ptrs = [1, 4, 7, 10, 13]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  end
  @notimplemented
end

"""
Provided the gids of the coarse cell faces as a Tuple(vertices,edges,cell) and the
edge-based RefinementRule we want to apply, returns the connectivity of all the refined
children cells.
"""
function get_relabeled_connectivity(rr::RefinementRule,faces_gids)
  return get_relabeled_connectivity(RefinementRuleType(rr),rr,faces_gids)
end

get_relabeled_connectivity(::RefinementRuleType,::RefinementRule,faces_gids) = @notimplemented

function get_relabeled_connectivity(::WithoutRefinement,rr::RefinementRule{P},faces_gids) where P<:Polytope{2}
  conn = rr.ref_grid.grid.cell_node_ids
  gids = faces_gids[1]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function get_relabeled_connectivity(::RedRefinement,rr::RefinementRule{P},faces_gids) where P<:Polytope{2}
  conn = rr.ref_grid.grid.cell_node_ids
  gids = [faces_gids[1]...,faces_gids[2]...,faces_gids[3]...]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function get_relabeled_connectivity(::GreenRefinement{N},rr::RefinementRule{P},faces_gids) where {N,P<:Polytope{2}}
  conn = rr.ref_grid.grid.cell_node_ids
  gids = [faces_gids[1]...,faces_gids[2][N]]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function get_relabeled_connectivity(::BlueRefinement{N, M},rr::RefinementRule{P},faces_gids) where {N, M, P<:Polytope{2}}
  conn = rr.ref_grid.grid.cell_node_ids
  gids = [faces_gids[1]...,faces_gids[2][N], faces_gids[2][M]]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function get_relabeled_connectivity(::BlueDoubleRefinement{N},rr::RefinementRule{P},faces_gids) where {N, P<:Polytope{2}}
  conn = rr.ref_grid.grid.cell_node_ids
  # Sort needed for non-longest edge gids
  other_ids = setdiff(collect(1:3), N) |> sort
  gids = [faces_gids[1]..., faces_gids[2][other_ids]..., faces_gids[2][N]]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function counts_to_ptrs(counts::Vector{<:Integer})
  n = length(counts)
  ptrs = Vector{Int32}(undef,n+1)
  @inbounds for i in 1:n
    ptrs[i+1] = counts[i]
  end
  length_to_ptrs!(ptrs)
  ptrs
end

