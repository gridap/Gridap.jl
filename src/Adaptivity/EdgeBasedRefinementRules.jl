
abstract type EdgeBasedRefinement <: RefinementRuleType end
struct RedRefinement <: EdgeBasedRefinement end
struct GreenRefinement{P} <: EdgeBasedRefinement end


"""
RefinementRule representing a non-refined cell.
"""
function WhiteRefinementRule(p::Polytope)
  ref_grid = UnstructuredDiscreteModel(UnstructuredGrid(LagrangianRefFE(Float64,p,1)))
  return RefinementRule(WithoutRefinement(),p,ref_grid)
end

function _has_interior_point(::RefinementRule{T}) where {T<:WithoutRefinement}
  return false
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
                4,2,6,
                5,6,3,
                6,5,4]
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

function _has_interior_point(rr::RefinementRule{T}) where {T<:RedRefinement}
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
  P = _get_green_vertex_permutation(p,ref_edge)
  if p == TRI
    polys     = [TRI]
    cell_type = [1,1]
    conn_data = [4,P[1],P[2],
                 4,P[3],P[1]]
    conn_ptrs = [1,4,7]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  elseif p == QUAD
    polys     = [TRI]
    cell_type = [1,1,1]
    conn_data = [P[1],P[2],5,
                 P[3],P[1],5,
                 P[4],5,P[2]]
    conn_ptrs = [1,4,7,10]
    return polys, cell_type, Table(conn_data,conn_ptrs)
  end
  @notimplemented
end

function _get_green_vertex_permutation(p::Polytope{2},ref_edge::Integer)
  if p == TRI
    perm = circshift([1,2,3],3-ref_edge)
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

function _has_interior_point(::RefinementRule{T}) where {T<:GreenRefinement}
  return false
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

function setup_edge_based_rrules(topo::UnstructuredGridTopology{Dc},cells_to_refine::Vector{<:Integer}) where Dc
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
      else # > 1, can't be 0
        cell_color[nbor] = RED + Int8(cell_types[nbor]-1)
        is_red[nbor] = true
        is_refined[nbor_edges] .= true
        enqueue!(q,nbor)
      end
    end
  end

  # Create RefinementRules
  num_rr = nP + nP + green_offsets[nP+1]-1 # WHITE+RED+GREEN
  color_rrules = Vector{RefinementRule}(undef,num_rr)
  for (k,p) in enumerate(polys)
    color_rrules[WHITE+k-1] = WhiteRefinementRule(p)
    color_rrules[RED+k-1]   = RedRefinementRule(p)
    for e in 1:num_faces(p,1)
      println(GREEN+green_offsets[k]-1+e-1)
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

  if all(lazy_map(rr->RefinementRuleType(rr)<:RedRefinement,rrules))
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

"""
Provided the gids of the coarse cell faces as a Tuple(vertices,edges,cell) and the 
edge-based RefinementRule we want to apply, returns the connectivity of all the refined 
children cells. 
"""
function get_relabeled_connectivity(rr::RefinementRule{T,P},faces_gids) where {P<:Polytope{2},T<:WithoutRefinement}
  conn = rr.ref_grid.grid.cell_node_ids
  gids = faces_gids[1]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function get_relabeled_connectivity(rr::RefinementRule{T,P},faces_gids) where {P<:Polytope{2},T<:RedRefinement}
  conn = rr.ref_grid.grid.cell_node_ids
  gids = [faces_gids[1]...,faces_gids[2]...,faces_gids[3]...]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function get_relabeled_connectivity(rr::RefinementRule{T,P},faces_gids) where {N,P<:Polytope{2},T<:GreenRefinement{N}}
  conn = rr.ref_grid.grid.cell_node_ids
  gids = [faces_gids[1]...,faces_gids[2][N]]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function counts_to_ptrs(counts)
  n = length(counts)
  ptrs = Vector{Int32}(undef,n+1)
  @inbounds for i in 1:n
    ptrs[i+1] = counts[i]
  end
  length_to_ptrs!(ptrs)
  ptrs
end

