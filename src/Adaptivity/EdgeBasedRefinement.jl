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

struct EdgeBasedRefinement <: AdaptivityMethod end

function refine(::EdgeBasedRefinement,model::UnstructuredDiscreteModel{Dc,Dp};cells_to_refine=nothing) where {Dc,Dp}
  # cells_to_refine can be 
  #    a) nothing -> All cells get refined
  #    b) AbstractArray{<:Bool} of size num_cells(model) 
  #            -> Only cells such that cells_to_refine[iC] == true get refined
  #    c) AbstractArray{<:Integer} 
  #            -> Cells for which gid ∈ cells_to_refine get refined
  ctopo = get_grid_topology(model)
  coarse_labels = get_face_labeling(model)

  # Create new model
  rrules, faces_list = setup_edge_based_rrules(ctopo,cells_to_refine)
  topo   = _refine_unstructured_topology(ctopo,rrules,faces_list)
  reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  grid   = UnstructuredGrid(get_vertex_coordinates(topo),get_faces(topo,Dc,0),reffes,get_cell_type(topo),OrientationStyle(topo))
  
  glue = _get_refinement_glue(topo,model.grid_topology,rrules)
  
  labels = _refine_face_labeling(coarse_labels,glue,model.grid_topology,topo)
  ref_model = UnstructuredDiscreteModel(grid,topo,labels)

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

The method returns a vector of Red/Green refinement rules, as well a the list of 
vertices, edges and cells which correspond to new vertices in the refined topology.
"""
function setup_edge_based_rrules(topo::UnstructuredGridTopology{Dc},::Nothing) where Dc
  rrules     = lazy_map(RedRefinementRule,CompressedArray(topo.polytopes,topo.cell_type))
  faces_list = _redgreen_refined_faces_list(topo,rrules,[true])
  return rrules, faces_list
end

function setup_edge_based_rrules(topo::UnstructuredGridTopology{Dc},cells_to_refine::AbstractArray{<:Bool}) where Dc
  return setup_edge_based_rrules(topo,findall(cells_to_refine))
end

function setup_edge_based_rrules(topo::UnstructuredGridTopology{Dc},cells_to_refine::AbstractArray{<:Integer}) where Dc
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
function RedRefinementRule(p::Polytope)
  @notimplementedif (p ∉ [TRI,QUAD,HEX])

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
  if p ∈ [QUAD,HEX]
    return true
  end
  return false
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
      Int32[1,2,3,1,2,3,1],           # [Fine node] -> [Coarse face]
      Int32[1,2,1,1,3,1,2,3,1],       # [Fine edge] -> [Coarse face]
      Int32[1,1,1,1]                  # [Fine cell] -> [Coarse face]
    ]
    parent_dims  = [
      Int32[0,0,0,1,1,1,2],           # [Fine node] -> [Coarse face dim]
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
    @notimplemented
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

function counts_to_ptrs(counts::Vector{<:Integer})
  n = length(counts)
  ptrs = Vector{Int32}(undef,n+1)
  @inbounds for i in 1:n
    ptrs[i+1] = counts[i]
  end
  length_to_ptrs!(ptrs)
  ptrs
end

