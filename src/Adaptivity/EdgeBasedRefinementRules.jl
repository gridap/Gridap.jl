
abstract type EdgeBasedRefinement <: RefinementRuleType end
struct RedRefinement <: EdgeBasedRefinement end
struct GreenRefinement{P} <: EdgeBasedRefinement end

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

  ref_grid = UnstructuredGrid(coords,conn,reffes,cell_types)
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

  ref_grid = UnstructuredGrid(coords,conn,reffes,cell_types)
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