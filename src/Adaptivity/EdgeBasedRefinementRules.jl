
abstract type EdgeBasedRefinement <: RefinementRuleType end
struct RedRefinement <: EdgeBasedRefinement end
struct GreenRefinement{P} <: EdgeBasedRefinement end

function WhiteRefinementRule(p::Polytope)
  ref_grid = UnstructuredDiscreteModel(UnstructuredGrid(LagrangianRefFE(Float64,p,1)))
  return RefinementRule(WithoutRefinement(),p,ref_grid)
end

function RedRefinementRule(p::Polytope{2})
  @notimplementedif (p ∉ [TRI,QUAD])

  if p == TRI
    faces_list = ([1,2,3],[1,2,3],[])
  elseif p == QUAD
    faces_list = ([1,2,3,4],[1,2,3,4],[1])
  end
  coords = _get_new_coordinates_from_faces(p,faces_list)

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

function GreenRefinementRule(p::Polytope{2},ref_edge::Integer)
  @notimplementedif (p ∉ [TRI,QUAD])

  faces_list = (collect(1:num_faces(p,0)),[ref_edge],[])
  coords = _get_new_coordinates_from_faces(p,faces_list)

  polys, cell_types, conn = _get_green_refined_connectivity(p,ref_edge)
  reffes = map(x->LagrangianRefFE(Float64,x,1),polys)

  ref_grid = UnstructuredDiscreteModel(UnstructuredGrid(coords,conn,reffes,cell_types))
  return RefinementRule(GreenRefinement{ref_edge}(),p,ref_grid)
end

function _get_green_refined_connectivity(p::Polytope{2},ref_edge)
  P = _get_green_permutation(p,ref_edge)
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

function _get_green_permutation(p::Polytope{2},ref_edge)
  e2n_map = get_faces(p,1,0)
  Ve = e2n_map[ref_edge]
  Vf = filter(x-> x ∉ Ve,collect(1:num_faces(p,0)))

  if p == TRI
    perm = [Vf[1],Ve[1],Ve[2]]
  elseif p == QUAD
    perm = [Vf[1],Vf[2],Ve[1],Ve[2]]
  else
    @notimplemented
  end

  return perm
end

function _get_new_coordinates_from_faces(p::Union{Polytope{D},GridTopology{D,Dp}},faces_list::Tuple) where {D,Dp}
  @check length(faces_list) == D+1

  nN_new     = sum(x->length(x),faces_list)
  coords_old = get_vertex_coordinates(p)
  coords_new = Vector{eltype(coords_old)}(undef,nN_new)

  n = 1
  for (d,dfaces) in enumerate(faces_list)
    if length(dfaces) > 0
      nf = length(dfaces)
      d2n_map = get_faces(p,d-1,0)
      nn = length(first(d2n_map))
      coords_new[n:n+nf-1] .= map(f -> sum(coords_old[d2n_map[f]])/nn, dfaces)
      n += nf
    end
  end

  return coords_new
end

function _setup_redgreen_coloring(topo::UnstructuredGridTopology{Dc},cells_to_refine::Vector{<:Integer}) where Dc
  nC = num_cells(topo)
  nE = num_faces(topo,1)
  c2e_map     = get_faces(topo,Dc,0)
  e2c_map     = get_faces(topo,0,Dc)
  polys       = topo.polytopes
  cell_types  = topo.cell_type
  nP          = length(polys)

  # Color pointers/offsets
  WHITE, RED, GREEN = Int8(1), Int8(1+nP), Int8(1+nP+nP)
  green_offsets = counts_to_ptrs(map(num_faces(p,1),polys))

  cell_color  = copy(cell_types) # WHITE
  is_red      = fill(false,nC)
  is_refined  = fill(false,nE)

  # Flag initial red cells and edges
  for c in cells_to_refine
    cell_color[c] .= RED + Int8(cell_types[c]-1)
    is_red[c] = true
    is_refined[c2e_map[c]] .= true
  end

  # Propagate red/green flags
  # Queue invariant: Cells in queue are RED, every RED cell is only touched once
  q = Queue{Int}()
  map!(c->enqueue!(q,c),cells_to_refine)
  while !isempty(q)
    c = dequeue!(q)
    c_edges = c2e_map[c]

    # For each non-red neighboring cell
    nbors = filter(n->!is_red[n],map(e->first(filter(x->x!=c,e2c_map[e])),c_edges))
    for nbor in nbors
      nbor_edges = c2e_map[nbor]
      nbor_ref_edges = findall(is_refined[nbor_edges])
      if length(nbor_ref_edges) == 1
        ref_edge = nbor_ref_edges[1]
        cell_color[nbor] = GREEN + Int8(green_offsets[ref_edge]+ref_edge-1)
      else # > 1, can't be 0
        cell_color[nbor] = RED + Int8(cell_types[nbor]-1)
        is_refined[nbor_edges] .= true
        enqueue!(q,nbor)
      end
    end
  end

  # Create RefinementRules
  num_rr = 1 + nP + green_offsets[nP+1]-1 # WHITE+RED+GREEN
  color_rrules = Vector{RefinementRule}(undef,num_rr)
  for (k,p) in enumerate(polys)
    color_rrules[WHITE+k-1] = WhiteRefinementRule(p)
    color_rrules[RED+k-1]   = RedRefinementRule(p)
    for e in 1:num_faces(p,1)
      color_rrules[GREEN+green_offsets[k]+e-1] = GreenRefinementRule(p,e)
    end
  end
  rrules = CompressedArray(color_rrules,cell_color)

  # Create list of refined faces
  faces_list = _redgreen_refined_faces_list(topo,rrules,is_refined)

  return rrules, faces_list
end

function _redgreen_refined_faces_list(topo::UnstructuredGridTopology{2},rrules,edge_is_refined)
  ref_nodes = 1:num_faces(topo,0)
  ref_edges = findall(edge_is_refined)
  ref_cells = findall(lazy_map(_has_interior_point,rrules))
  return (ref_nodes,ref_edges,ref_cells)
end

function _has_interior_point(rr::RefinementRule{P,T}) where {P,T<:RedRefinement}
  p = get_polytope(rr)
  if p == QUAD
    return true
  end
  return false
end

function _has_interior_point(rr::RefinementRule{P,T}) where {P,T<:GreenRefinement}
  return false
end

function get_relabeled_connectivity(rr::RefinementRule{P,T},faces_gids) where {P<:Polytope{2},T<:RedRefinement}
  conn = rr.ref_grid.grid.cell_node_ids
  gids = [faces_gids[1]...,faces_gids[2]...,faces_gids[3]...]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

function get_relabeled_connectivity(rr::RefinementRule{P,T},faces_gids) where {N,P<:Polytope{2},T<:GreenRefinement{N}}
  conn = rr.ref_grid.grid.cell_node_ids
  gids = [faces_gids[1]...,faces_gids[2][N]]
  new_data = lazy_map(Reindex(gids),conn.data)
  return Table(new_data,conn.ptrs)
end

"""
function _get_green_refined_connectivity_bis(poly::Polytope{2},ref_edge)
  nChildren = num_faces(poly,0)-1
  nData = 3*nChildren
  polys = Fill(TRI,nChildren)

  e2n_map = get_faces(poly,1,0)
  E_ref   = e2n_map[ref_edge]
  remaining_edges = filter(x->!(x==E_ref),e2n_map)

  conn_data = Vector{Int32}(undef,nData)
  conn_ptrs = Vector{Int32}(undef,nChildren+1)

  conn_ptrs[1] = 1
  for (k,E) in enumerate(remaining_edges)
    pk = conn_ptrs[k]
    conn_ptrs[k+1] = pk+3
    conn_data[pk:pk+2] .= [4,E[1],E[2]]
  end

  return Table(conn_data,conn_ptrs)
end
"""