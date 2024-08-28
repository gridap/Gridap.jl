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


"""
    abstract type EdgeBasedRefinement <: AdaptivityMethod end

Abstract type representing edge-based refinement methods. Edge-based refinement methods
are the ones where all new nodes are created from the faces of the parent mesh.
This allows the new connectivity to easily be computed from the parent connectivity.

Any edge-based refinement method has to implement the following API: 

- [`setup_edge_based_rrules`](@ref)

"""
abstract type EdgeBasedRefinement <: AdaptivityMethod end

function refine(method::EdgeBasedRefinement,model::UnstructuredDiscreteModel{Dc,Dp};cells_to_refine=nothing) where {Dc,Dp}
  # cells_to_refine can be
  #    a) nothing -> All cells get refined
  #    b) AbstractArray{<:Bool} of size num_cells(model)
  #            -> Only cells such that cells_to_refine[iC] == true get refined
  #    c) AbstractArray{<:Integer} 
  #            -> Cells for which gid ∈ cells_to_refine get refined
  ctopo = get_grid_topology(model)
  coarse_labels = get_face_labeling(model)
  # Create new model
  rrules, faces_list = setup_edge_based_rrules(method, model.grid_topology,cells_to_refine)
  topo   = refine_edge_based_topology(ctopo,rrules,faces_list)
  reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  grid   = UnstructuredGrid(
    get_vertex_coordinates(topo),
    get_faces(topo,Dc,0),reffes,
    get_cell_type(topo),
    OrientationStyle(topo)
  )
  glue = blocked_refinement_glue(rrules)
  labels = refine_face_labeling(coarse_labels,glue,model.grid_topology,topo)
  ref_model = UnstructuredDiscreteModel(grid,topo,labels)
  return AdaptedDiscreteModel(ref_model,model,glue)
end

function refine_edge_based_topology(
  topo::UnstructuredGridTopology{Dc},
  rrules::AbstractVector{<:RefinementRule},
  faces_list::Tuple
) where {Dc}
  coords_new  = get_new_coordinates_from_faces(topo,faces_list)
  c2n_map_new = get_refined_cell_to_vertex_map(topo,rrules,faces_list)
  polys_new, cell_type_new = _get_cell_polytopes(rrules)
  orientation = NonOriented()

  return UnstructuredGridTopology(coords_new,c2n_map_new,cell_type_new,polys_new,orientation)
end

"""
    setup_edge_based_rrules(method::EdgeBasedRefinement,topo::UnstructuredGridTopology,cells_to_refine)

Given an UnstructuredTopology and a list of cells to refine, returns a vector
of refinement rules and a list of the faces (called refined faces) where new 
vertices will be created.
"""
function setup_edge_based_rrules(::EdgeBasedRefinement,topo::UnstructuredGridTopology{Dc},cells_to_refine) where Dc
  @abstractmethod
end

function setup_edge_based_rrules(method::EdgeBasedRefinement,topo::UnstructuredGridTopology{Dc},cells_to_refine::AbstractArray{<:Bool}) where Dc
  return setup_edge_based_rrules(method,topo,findall(cells_to_refine))
end

"""
    get_new_coordinates_from_faces(p::Union{Polytope{D},GridTopology{D}},faces_list::Tuple)

Given a Polytope or GridTopology and a list of refined d-faces, returns the new coordinates
of the vertices created on those faces: 
  - If the face is a 0-face (node), the new vertex is the node itself.
  - For d > 0, the new vertex is the barycenter of the face.
The new vertices are ordered by parent dimension and face id (in that order).
"""
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

"""
    get_refined_cell_to_vertex_map(
      topo::UnstructuredGridTopology{Dc},
      rrules::AbstractVector{<:RefinementRule},
      faces_list::Tuple
    )

Given an UnstructuredGridTopology, a vector of RefinementRules and a list of refined faces,
returns the cell-to-vertex connectivity of the refined topology.
"""
function get_refined_cell_to_vertex_map(
  topo::UnstructuredGridTopology{Dc},
  rrules::AbstractVector{<:RefinementRule},
  faces_list::Tuple
) where Dc
  @check length(faces_list) == Dc+1
  d_to_nfaces = map(Df -> num_faces(topo,Df), 0:Dc)
  d_to_cell_to_dface = map(Df -> get_faces(topo,Dc,Df), 0:Dc-1)

  # Offsets for new vertices
  d_to_offset = [1,map(length,faces_list)...]
  length_to_ptrs!(d_to_offset)
  d_to_offset = d_to_offset[1:end-1]

  # Maps from refined faces (i.e new vertices) to old faces
  d_to_faces_reindexing = map(d_to_nfaces,d_to_offset,faces_list) do nfaces,offset,ref_faces
    faces_reindexing = find_inverse_index_map(ref_faces,nfaces)
    return map(f -> iszero(f) ? 0 : f + offset - 1, faces_reindexing)
  end

  # For a given parent cell, return the new node gids (corresponding to the refined faces)
  function cell_to_reindexed_faces(::Val{2},cell)
    N = view(d_to_cell_to_dface[1],cell)
    E = view(d_to_faces_reindexing[2],view(d_to_cell_to_dface[2],cell))
    C = [d_to_faces_reindexing[3][cell]]
    return (N,E,C)
  end
  function cell_to_reindexed_faces(::Val{3},cell)
    N = view(d_to_cell_to_dface[1],cell)
    E = view(d_to_faces_reindexing[2],view(d_to_cell_to_dface[2],cell))
    F = view(d_to_faces_reindexing[3],view(d_to_cell_to_dface[3],cell))
    C = [d_to_faces_reindexing[4][cell]]
    return (N,E,F,C)
  end

  # Allocate ptr and data arrays for new connectivity
  nC_new    = sum(rr -> num_subcells(rr), rrules)
  nData_new = sum(rr -> sum(length,rr.ref_grid.grid.cell_node_ids), rrules)
  ptrs_new  = Vector{Int}(undef,nC_new+1)
  data_new  = Vector{Int}(undef,nData_new)

  k = 1
  ptrs_new[1] = 1
  for iC = 1:d_to_nfaces[Dc+1]
    rr = rrules[iC]

    # New vertex ids from old face ids
    faces = cell_to_reindexed_faces(Val(Dc),iC)
    sub_conn = get_relabeled_connectivity(rr,faces)

    nChild = length(sub_conn)
    ptrs_new[k:k+nChild] .= sub_conn.ptrs .+ (ptrs_new[k] - 1)
    data_new[ptrs_new[k]:ptrs_new[k+nChild]-1] .= sub_conn.data
    k = k+nChild
  end

  return Table(data_new,ptrs_new)
end

############################################################################################

# NVB Refinement (Nearest Vertex Bissection)

"""
    struct NVBRefinement{T} <: EdgeBasedRefinement

A refinement method for TRIs, where a selected cell is refined by bisecting the longest edge.
"""
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

function setup_edge_based_rrules(method::NVBRefinement, topo::UnstructuredGridTopology{Dc},::Nothing) where Dc
  setup_edge_based_rrules(method, topo, collect(1:num_faces(topo,Dc)))
end

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

############################################################################################

# RedGreen Refinement

"""
    struct RedGreenRefinement <: EdgeBasedRefinement

Red-Green conformal refinement method. The refinement follows the following rules:

  - If a cell is completely refined, it is colored RED
  - If a cell touches more than one refined (RED) cell, it becomes RED.
  - If a cell touches a single refined (RED) cell, it is colored GREEN.

The resulting mesh is always conformal.
"""
struct RedGreenRefinement <: EdgeBasedRefinement end

function setup_edge_based_rrules(::RedGreenRefinement, topo::UnstructuredGridTopology{Dc},::Nothing) where Dc
  rrules = lazy_map(RedRefinementRule,CompressedArray(topo.polytopes,topo.cell_type))

  ref_nodes = 1:num_faces(topo,0)
  ref_edges = 1:num_faces(topo,1)
  ref_cells = findall(_has_interior_point,rrules)

  if Dc == 2
    faces_list = (ref_nodes,ref_edges,ref_cells)
  else
    @check length(topo.polytopes) == 1 "Mixed meshes not supported yet."
    p = first(topo.polytopes)
    ref_faces = is_simplex(p) ? Int32[] : 1:num_faces(topo,2)
    faces_list = (ref_nodes,ref_edges,ref_faces,ref_cells)
  end

  return rrules, faces_list
end

function setup_edge_based_rrules(
  ::RedGreenRefinement, topo::UnstructuredGridTopology{Dc},cells_to_refine::AbstractArray{<:Integer}
) where Dc
  @check Dc == 2
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

############################################################################################

# Barycentric refinement

"""
    struct BarycentricRefinementMethod <: EdgeBasedRefinement end

Barycentric refinement method. Each selected cell is refined by adding a new vertex at the
barycenter of the cell. Then the cell is split into subcells by connecting the new vertex
to the vertices of the parent cell. Always results in a conformal mesh.
"""
struct BarycentricRefinement <: EdgeBasedRefinement end

function setup_edge_based_rrules(
  method::BarycentricRefinement, topo::UnstructuredGridTopology{Dc},::Nothing
) where Dc
  setup_edge_based_rrules(method,topo,collect(1:num_faces(topo,Dc)))
end

function setup_edge_based_rrules(
  ::BarycentricRefinement, topo::UnstructuredGridTopology{Dc},cells_to_refine::AbstractArray{<:Integer}
) where Dc
  @check (length(topo.polytopes) == 1) && is_simplex(first(topo.polytopes)) "Only simplex meshes supported"

  ptrs = fill(1,num_faces(topo,Dc))
  ptrs[cells_to_refine] .= 2
  p = first(topo.polytopes)
  rrules = CompressedArray([WhiteRefinementRule(p),BarycentricRefinementRule(p)],ptrs)

  ref_nodes = 1:num_faces(topo,0)
  if Dc == 2
    faces_list = (ref_nodes,Int32[],cells_to_refine)
  else
    faces_list = (ref_nodes,Int32[],Int32[],cells_to_refine)
  end

  return rrules, faces_list
end

############################################################################################

# EdgeBasedRefinementRules

abstract type EdgeBasedRefinementRule <: RefinementRuleType end

# List of available RefinementRules:
struct RedRefinement <: EdgeBasedRefinementRule end
struct GreenRefinement{P} <: EdgeBasedRefinementRule end
struct BlueRefinement{P, Q} <: EdgeBasedRefinementRule end
struct BlueDoubleRefinement{P} <: EdgeBasedRefinementRule end
struct BarycentricRefinementRule <: EdgeBasedRefinementRule end
struct PowellSabinRefinement <: EdgeBasedRefinementRule end

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
function RedRefinementRule(p::Polytope)
  @notimplementedif (p ∉ [TRI,TET,QUAD,HEX])

  faces_list = _get_red_refined_faces_list(p)
  coords = get_new_coordinates_from_faces(p,faces_list)

  polys, cell_types, conn = _get_red_refined_connectivity(p)
  reffes = map(x->LagrangianRefFE(Float64,x,1),polys)

  ref_grid = UnstructuredDiscreteModel(UnstructuredGrid(coords,conn,reffes,cell_types))
  return RefinementRule(RedRefinement(),p,ref_grid)
end

function _get_red_refined_faces_list(p::Polytope)
  if p == TRI
    return (Int32[1,2,3],Int32[1,2,3],Int32[])
  elseif p == TET
    return (Int32[1,2,3,4],Int32[1,2,3,4,5,6],Int32[],Int32[])
  elseif p == QUAD
    return (Int32[1,2,3,4],Int32[1,2,3,4],Int32[1])
  elseif p == HEX
    return (Int32[1,2,3,4,5,6,7,8],Int32[1,2,3,4,5,6,7,8,9,10,11,12],Int32[1,2,3,4,5,6],Int32[1])
  end
  @notimplemented
end

function _get_red_refined_connectivity(p::Polytope)
  if p == TRI
    polys     = [TRI]
    cell_type = Int32[1,1,1,1]
    conn_data = Int32[1,4,5,
                      2,4,6,
                      3,5,6,
                      4,5,6]
    conn_ptrs = Int32[1,4,7,10,13]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  elseif p == TET
    polys     = [TET]
    cell_type = Int32[1,1,1,1,1,1,1,1]
    conn_data = Int32[1,5,6,8,
                      5,2,7,9,
                      6,7,3,10,
                      8,9,10,4,
                      5,7,10,9,
                      8,5,10,9,
                      5,10,7,6,
                      8,10,5,6]
    conn_ptrs = Int32[1,5,9,13,17,21,25,29,33]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  elseif p == QUAD
    polys     = [QUAD]
    cell_type = Int32[1,1,1,1]
    conn_data = Int32[1,5,7,9,
                      5,2,9,8,
                      7,9,3,6,
                      9,8,6,4]
    conn_ptrs = Int32[1,5,9,13,17]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  elseif p == HEX
    polys     = [HEX]
    cell_type = Int32[1,1,1,1,1,1,1,1]
    conn_data = Int32[ 1, 9,13,21,17,23,25,27,
                       9, 2,21,14,23,18,27,26,
                      13,21, 3,10,25,27,19,24,
                      21,14,10, 4,27,26,24,20,
                      17,23,25,27, 5,11,15,22,
                      23,18,27,26,11, 6,22,16,
                      25,27,19,24,15,22, 7,12,
                      27,26,24,20,22,16,12, 8]
    conn_ptrs = Int32[1,9,17,25,33,41,49,57,65]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  end
  @notimplemented
end

function _has_interior_point(rr::RefinementRule,::RedRefinement)
  p = get_polytope(rr)
  return !is_simplex(p)
end

# [Face dimension][Coarse Face id] -> [Fine faces]
function get_d_to_face_to_child_faces(rr::RefinementRule,::RedRefinement)
  p = get_polytope(rr)
  if p == QUAD
    return [
      [Int32[1],Int32[2],Int32[3],Int32[4]],           # [Coarse Node] -> [Fine Node]
      [Int32[1,5],Int32[8,11],Int32[3,9],Int32[7,12]], # [Coarse Edge] -> [Fine Edge]
      [Int32[1,2,3,4]]                                 # [Coarse Cell] -> [Fine Cells]
    ]
  elseif p == TRI
    return [
      [Int32[1],Int32[2],Int32[3]],       # [Coarse Node] -> [Fine Node]
      [Int32[1,4],Int32[2,7],Int32[5,8]], # [Coarse Edge] -> [Fine Edge]
      [Int32[1,2,3]]                      # [Coarse Cell] -> [Fine Cells]
    ]
  elseif p == HEX
    return [
      [Int32[1],Int32[2],Int32[3],Int32[4],
       Int32[5],Int32[6],Int32[7],Int32[8]],                # [Coarse Node] -> [Fine Node]
      [Int32[1,13],Int32[21,29],Int32[34,42],Int32[47,52],
       Int32[5,23],Int32[17,31],Int32[36,48],Int32[44,53],
       Int32[9,38],Int32[19,45],Int32[27,50],Int32[33,54]], # [Coarse Edge] -> [Fine Edge]
      [Int32[1,7,12,17] ,Int32[21,26,30,34],
       Int32[3,9,22,27] ,Int32[14,19,31,35],
       Int32[5,15,24,32],Int32[11,20,29,36]],               # [Coarse Face] -> [Fine Face]
      [Int32[1,2,3,4,5,6,7,8]],                             # [Coarse Cell] -> [Fine Cells]
    ]
  else
    @notimplemented
  end
end

# 1 - [Face dimension][Fine Face id] -> [Parent Face]
# 2 - [Face dimension][Fine Face id] -> [Parent Face Dimension]
function get_d_to_face_to_parent_face(rr::RefinementRule,::RedRefinement)
  p = get_polytope(rr)
  if p == QUAD
    parent_faces = [
      Int32[1,2,3,4,1,2,3,4,1],       # [Fine node] -> [Coarse face]
      Int32[1,1,3,1,1,1,4,2,3,1,2,4], # [Fine edge] -> [Coarse face]
      Int32[1,1,1,1]                  # [Fine cell] -> [Coarse face]
    ]
    parent_dims  = [
      Int32[0,0,0,0,1,1,1,1,2],       # [Fine node] -> [Coarse face dim]
      Int32[1,2,1,2,1,2,1,1,1,2,1,1], # [Fine edge] -> [Coarse face dim]
      Int32[2,2,2,2]                  # [Fine cell] -> [Coarse face dim]
    ]
  elseif p == TRI
    parent_faces = [
      Int32[1,2,3,1,2,3],             # [Fine node] -> [Coarse face]
      Int32[1,2,1,1,3,1,2,3,1],       # [Fine edge] -> [Coarse face]
      Int32[1,1,1,1]                  # [Fine cell] -> [Coarse face]
    ]
    parent_dims  = [
      Int32[0,0,0,1,1,1],             # [Fine node] -> [Coarse face dim]
      Int32[1,1,2,1,1,2,1,1,2],       # [Fine edge] -> [Coarse face dim]
      Int32[2,2,2,2]                  # [Fine cell] -> [Coarse face dim]
    ]
  elseif p == HEX
    parent_faces = [ 
      Int32[1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,1], # [Fine node] -> [Coarse face]
      Int32[1,1,3,1,5,1,5,1,9,3,5,1,
            1,1,3,1,6,6,10,6,
            2,4,5,1,5,1,11,4,
            2,4,6,6,12,
            3,2,7,2,9,3,5,1,
            3,2,8,10,6,
            4,7,2,11,4,
            4,8,12],         # [Fine edge] -> [Coarse face]
      Int32[1,1,3,1,5,1,
            1,1,3,1,6,
            1,1,4,5,1,
            1,1,4,6,
            2,3,1,5,1,
            2,3,1,6,
            2,4,5,1,
            2,4,6],          # [Fine face] -> [Coarse face]
      Int32[1,1,1,1,1,1,1,1] # [Fine cell] -> [Coarse face]
    ]
    parent_dims = [
      Int32[0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3], # [Fine node] -> [Coarse face dim]
      Int32[1,2,2,3,1,2,2,3,1,2,2,3,
            1,2,2,3,1,2,1,2,
            1,2,1,2,2,3,1,2,
            1,2,1,2,1,
            1,2,1,2,1,2,2,3,
            1,2,1,1,2,
            1,1,2,1,2,
            1,1,1],          # [Fine edge] -> [Coarse face dim]
      Int32[2,3,2,3,2,3,
            2,3,2,3,2,
            2,3,2,2,3,
            2,3,2,2,
            2,2,3,2,3,
            2,2,3,2,
            2,2,2,3,
            2,2,2],          # [Fine face] -> [Coarse face dim]
      Int32[3,3,3,3,3,3,3,3] # [Fine cell] -> [Coarse face dim]
    ]
  else
    parent_faces, parent_dims = get_d_to_face_to_parent_face(rr,GenericRefinement())
  end
  return parent_faces, parent_dims
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
RefinementRule representing barycentric refinement of a Polytope. A single
vertex is added in the center, then joined to the vertices of the original
Polytope.
"""
function BarycentricRefinementRule(p::Polytope)
  @notimplementedif (p ∉ [TRI,TET])

  faces_list = _get_barycentric_refined_faces_list(p)
  coords = get_new_coordinates_from_faces(p,faces_list)

  polys, cell_types, conn = _get_barycentric_refined_connectivity(p)
  reffes = map(x->LagrangianRefFE(Float64,x,1),polys)

  ref_grid = UnstructuredDiscreteModel(UnstructuredGrid(coords,conn,reffes,cell_types))
  return RefinementRule(BarycentricRefinementRule(),p,ref_grid)
end

function _get_barycentric_refined_faces_list(p::Polytope)
  if p == TRI
    return (Int32[1,2,3],Int32[],Int32[1])
  elseif p == TET
    return (Int32[1,2,3,4],Int32[],Int32[],Int32[1])
  end
  @notimplemented
end

function _get_barycentric_refined_connectivity(p::Polytope)
  if p == TRI
    polys     = [TRI]
    cell_type = [1, 1, 1]
    conn_data = [1, 2, 4,
                 2, 3, 4,
                 3, 1, 4]
    conn_ptrs = [1, 4, 7, 10]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  elseif p == TET
    polys     = [TET]
    cell_type = [1, 1, 1, 1]
    conn_data = [1, 2, 3, 5,
                 2, 3, 4, 5,
                 3, 1, 4, 5,
                 1, 2, 4, 5]
    conn_ptrs = [1, 5, 9, 13, 17]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  end
  @notimplemented
end

"""
RefinementRule representing the Powell-Sabin split of TRI, and it's extension to
3-demensional TETs (Worsey-Farin split).
"""
function PowellSabinRefinementRule(p::Polytope)
  @notimplementedif (p ∉ [TRI,TET])

  faces_list = _get_powellsabin_refined_faces_list(p)
  coords = get_new_coordinates_from_faces(p,faces_list)

  polys, cell_types, conn = _get_powellsabin_refined_connectivity(p)
  reffes = map(x->LagrangianRefFE(Float64,x,1),polys)

  ref_grid = UnstructuredDiscreteModel(UnstructuredGrid(coords,conn,reffes,cell_types))
  return RefinementRule(PowellSabinRefinement(),p,ref_grid)
end

function _get_powellsabin_refined_faces_list(p::Polytope)
  if p == TRI
    return (Int32[1,2,3],Int32[1,2,3],Int32[1])
  elseif p == TET
    return (Int32[1,2,3,4],Int32[1,2,3,4,5,6],Int32[1,2,3,4],Int32[1])
  end
  @notimplemented
end

function _get_powellsabin_refined_connectivity(p::Polytope)
  # A) For TRIs, we add a new vertices 
  #       - in the center of the triangle
  #       - in the middle of each edge
  # and then connect the center vertex to all other vertices (nodes + new edge vertices). 
  # This yields 6 new triangles.
  # B) For TETs, we do do the same process for each face (which is a TRI), and then connect 
  #    every new triangle to the center of the TET. This yields 24=n_TRI*n_FACE new TETs.

  n_TRI  = 6 # Number of children in a refined TRI
  n_FACE = 4 # Number of faces in a TET
  conn_data_TRI = [1, 4, 7,
                  4, 2, 7,
                  2, 6, 7,
                  6, 3, 7,
                  3, 5, 7,
                  5, 1, 7]
  if p == TRI
    polys = [TRI]
    cell_type = fill(1, n_TRI)
    conn_ptrs = collect(1:3:n_TRI*3+1)
    return polys, cell_type, Table(conn_data_TRI,conn_ptrs)
  elseif p == TET
    polys     = [TET]
    cell_type = fill(1, n_TRI*n_FACE)
    face_conn = [[1,2,3,5,6,7,11],[1,2,4,5,8,9,12],[1,3,4,6,8,10,13],[2,3,4,7,9,10,14]]
    conn_data = vcat([[lazy_map(Reindex(face_conn[f]),conn_data_2d[1+(i-1)*3:i*3])...,15] for i in 1:n_TRI for f in 1:n_FACE]...)
    conn_ptrs = collect(1:4:n_TRI*n_FACE*4+1)
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

function get_relabeled_connectivity(::WithoutRefinement,rr::RefinementRule,faces_gids)
  conn = rr.ref_grid.grid.cell_node_ids
  gids = faces_gids[1]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function get_relabeled_connectivity(::RedRefinement,rr::RefinementRule,faces_gids)
  conn = rr.ref_grid.grid.cell_node_ids
  if is_simplex(get_polytope(rr)) # TRI, TET only adds vertices to edges
    gids = vcat(faces_gids[1]...,faces_gids[2]...)
  else # QUAD, HEX adds vertices to edges, faces and interior
    gids = vcat(faces_gids...)    
  end
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function get_relabeled_connectivity(::GreenRefinement{N},rr::RefinementRule{<:Polytope{2}},faces_gids) where N
  conn = rr.ref_grid.grid.cell_node_ids
  gids = [faces_gids[1]...,faces_gids[2][N]]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function get_relabeled_connectivity(::BlueRefinement{N, M},rr::RefinementRule{<:Polytope{2}},faces_gids) where {N, M}
  conn = rr.ref_grid.grid.cell_node_ids
  gids = [faces_gids[1]...,faces_gids[2][N], faces_gids[2][M]]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function get_relabeled_connectivity(::BlueDoubleRefinement{N},rr::RefinementRule{<:Polytope{2}},faces_gids) where {N}
  conn = rr.ref_grid.grid.cell_node_ids
  # Sort needed for non-longest edge gids
  other_ids = setdiff(collect(1:3), N) |> sort
  gids = [faces_gids[1]..., faces_gids[2][other_ids]..., faces_gids[2][N]]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function get_relabeled_connectivity(::BarycentricRefinementRule,rr::RefinementRule{<:Polytope{D}},faces_gids) where D
  conn = rr.ref_grid.grid.cell_node_ids
  gids = [faces_gids[1]...,faces_gids[D+1]...] # Polytope nodes + Barycenter
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function get_relabeled_connectivity(::PowellSabinRefinement,rr::RefinementRule,faces_gids)
  @check is_simplex(get_polytope(rr))
  conn = rr.ref_grid.grid.cell_node_ids
  gids = vcat(faces_gids...)
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

