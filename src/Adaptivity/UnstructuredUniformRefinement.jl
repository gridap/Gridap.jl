@inline function _point_isless(x::Point{D},y::Point{D},atol=eps()) where {D}
  @inbounds for i ∈ D:-1:1
    if x[i] < y[i] - atol
      return true
    elseif x[i] > y[i] + atol
      return false
    end
  end
  true
end

@inline function _to_cube_pattern(x::Point,atol=eps())
  map(x.data) do xi
    if isapprox(xi,zero(xi);atol=atol)
      -one(Int8)
    elseif isapprox(xi,one(xi),atol=atol)
      zero(Int8)
    else
      one(Int8)
    end
  end
end

@inline function _to_simplex_pattern(x::Point,atol=eps())
  flag = isapprox(sum(x),one(eltype(x));atol=atol)
  (_to_cube_pattern(x,atol)..., flag ? one(Int8) : zero(Int8))
end

function _num_interior_nodes(p::Polytope,n)::Int
  Dc = num_cell_dims(p)
  if is_n_cube(p)
    _num_nodes(p,n-2)
  elseif is_simplex(p)
    _num_nodes(p,n-Dc-1)
  else
    @notimplemented
  end
end

function _num_cells(p::Polytope,n)::Int
  @notimplementedif !(is_n_cube(p) || is_simplex(p))
  n^(num_dims(p))
end

function _num_nodes(p::Polytope,n)::Int
  Dc = num_cell_dims(p)
  if is_n_cube(p)
    (n+1)^Dc
  elseif is_simplex(p)
    binomial(Dc+n,n)
  else
    @notimplemented
  end
end

function cube_simplex_pattern_dimfid(p::Polytope{D}) where {D}
  face_bcs::Vector{Point{D}} = map(mean,get_face_coordinates(p))
  dimfid = mapreduce(vcat,0:D) do d
    broadcast(tuple,d,eachindex(get_dimrange(p,d)))
  end
  if is_n_cube(p)
    pt = map(_to_cube_pattern,face_bcs)
  else
    pt = map(_to_simplex_pattern,face_bcs)
  end
  TKey = Union{NTuple{D,Int8},NTuple{D+1,Int8}}
  Dict{TKey,Tuple{Int,Int}}(broadcast(Pair,pt,dimfid)...)
end

function cube_simplex_reference_grid(p::Polytope,n::Integer,pt_map)
  Dc = num_cell_dims(p)
  atol = 0.25 / n 
  if is_n_cube(p)
    ref_grid = UnstructuredGrid(Geometry._compute_linear_grid_from_n_cube(p,n))
    coords = get_node_coordinates(ref_grid)
    dimfids = map(c->pt_map[_to_cube_pattern(c,atol)],coords)
  else 
    ref_grid = Geometry._compute_linear_grid_from_simplex(p,n)
    coords = get_node_coordinates(ref_grid)
    dimfids = map(c->pt_map[_to_simplex_pattern(c,atol)],coords)
  end
  
  perm = sortperm(broadcast(tuple,dimfids,coords))
  conn = Table(get_cell_node_ids(ref_grid))
  inv_perm = similar(perm)
  inv_perm[perm] = 1:length(perm)
  new_coords = coords[perm]
  new_conn = Table(inv_perm[conn.data],conn.ptrs)
  reffes::Vector{LagrangianRefFE{Dc}} = get_reffes(ref_grid)
  orientation_style = NonOriented()
  cell_types = collect(get_cell_type(ref_grid))

  UnstructuredGrid(
    new_coords,
    new_conn,
    reffes,
    cell_types,
    orientation_style
  )
end

function compute_d_dface_offsets(ctopo,n,cell_refine_masks::AbstractVector{Bool})
  @assert length(cell_refine_masks) == num_cells(ctopo)
  
  Dc = num_cell_dims(ctopo)
  cell_dface_nnodes = CompressedArray(
    map(
      Broadcasting(p->_num_interior_nodes(p,n)),
      map(get_reffaces,get_polytopes(ctopo))
    )::Vector{Vector{Int}},
    get_cell_type(ctopo)
  )
  d_dface_offsets = Vector{Vector{Int32}}(undef,Dc)
  for d in 1:Dc
    (;data,ptrs) = Table(get_faces(ctopo,d,Dc))
    d_offsets = similar(ptrs,Int32)
    d_offsets[1] = d == 1 ? num_vertices(ctopo)+1 : d_dface_offsets[d-1][end]
    @inbounds for dfi = 1:length(d_offsets)-1
      I = findfirst(i->cell_refine_masks[i],view(data,ptrs[dfi]:ptrs[dfi+1]-1))
      if !isnothing(I)
        ci = data[ptrs[dfi]] + I - 1
        d_offsets[dfi+1] = d_offsets[dfi] + cell_dface_nnodes[ci][d+1]
      else
        d_offsets[dfi+1] = d_offsets[dfi]
      end
    end
    d_dface_offsets[d] = d_offsets
  end
  d_dface_offsets
end

function compute_cell_offsets(ctopo,n,cell_refine_masks::AbstractVector{Bool})
  @assert length(cell_refine_masks) == num_cells(ctopo)

  Dc = num_cell_dims(ctopo)
  polytopes = get_polytopes(ctopo)
  ptrs = collect(get_cell_type(ctopo))
  values = vcat(
    map(num_vertices,polytopes),
    map(p->_num_nodes(p,n),polytopes)
  )
  ptrs[cell_refine_masks] .+= length(polytopes)
  cell_nnodes = CompressedArray(values,ptrs)
  offsets = Vector{Int32}(undef,num_faces(ctopo,Dc)+1)
  offsets[1] = 1
  @inbounds for ci in 1:length(offsets)-1
    offsets[ci+1] = offsets[ci] + cell_nnodes[ci]
  end
  offsets
end

function unstructured_uniform_cell_l2gmap_and_nnodes(ctopo,n,cell_refine_masks::AbstractVector{Bool})
  @assert length(cell_refine_masks) == num_cells(ctopo)

  Dc = num_cell_dims(ctopo)
  d_df_goffsets = compute_d_dface_offsets(ctopo,n,cell_refine_masks)
  l2g_ptrs = compute_cell_offsets(ctopo,n,cell_refine_masks)
  l2g_data = Vector{Int32}(undef,l2g_ptrs[end]-1)
  d_c2df = ntuple(d->Table(get_faces(ctopo,Dc,d)),Val{Dc}())
  c2n = Table(get_faces(ctopo,Dc,0))

  @inbounds Threads.@threads for ci in 1:length(l2g_ptrs)-1
    lo = l2g_ptrs[ci]

    # vertices
    (;data,ptrs) = c2n
    l = ptrs[ci+1]-ptrs[ci]
    copyto!(l2g_data,lo,data,ptrs[ci],l)
    lo += l

    # 1-face ~ (Dc)-face
    if cell_refine_masks[ci]
      for d in 1:Dc
        (;data,ptrs) = d_c2df[d]
        goffsets = d_df_goffsets[d]
        for i in ptrs[ci]:ptrs[ci+1]-1
          fi = data[i]
          for go in goffsets[fi]:goffsets[fi+1]-1
            l2g_data[lo] = go
            lo += 1
          end
        end
      end
    end
  end

  n_nodes = d_df_goffsets[end][end] - 1
  Table(l2g_data,l2g_ptrs),n_nodes
end

function unstructured_uniform_coordinates(cell_ref_coords,cell_l2g,cell_map,n_nodes)
  coords = similar(first(cell_ref_coords),n_nodes)
  cell_l2g = Table(cell_l2g)
  @inbounds Threads.@threads for ci in eachindex(cell_l2g)
    ref_co = cell_ref_coords[ci]
    (;data,ptrs) = cell_l2g
    f = cell_map[ci]
    pini = ptrs[ci]
    l = ptrs[ci+1]-pini
    for j in 1:l
      coords[data[pini+j-1]] = f(ref_co[j])
    end
  end
  coords
end

function unstructured_uniform_connectivity(cell_ref_conns,cell_l2g,n_cells)
  cell_ref_conns = lazy_map(Table, cell_ref_conns)
  conn_ptrs = Vector{Int32}(undef,n_cells+1)
  conn_ptrs[1] = 1
  fci = 1
  @inbounds for ci in eachindex(cell_ref_conns)
    (;ptrs) = cell_ref_conns[ci]
    # cptrs = conn.ptrs
    for i in 1:length(ptrs)-1
      conn_ptrs[fci+1] = conn_ptrs[fci] + ptrs[i+1] - ptrs[i]
      fci += 1
    end
  end

  conn_data = Vector{Int32}(undef,conn_ptrs[end]-1)
  cache = array_cache(cell_l2g)
  go = 1
  @inbounds for ci in eachindex(cell_ref_conns)
    (;data) = cell_ref_conns[ci]
    n = length(data)
    l2g = getindex!(cache,cell_l2g,ci)
    conn_data[go:go+n-1] = l2g[data]
    go += n
  end
  Table(conn_data,conn_ptrs)
end

function unstructured_uniform_topology(cell_ref_grid,ctopo,cell_map,n,cell_refine_masks)
  polytopes = get_polytopes(ctopo)
  @notimplementedif !all( map(p->is_n_cube(p) || is_simplex(p),polytopes) )

  cell_lncells = lazy_map(num_cells,cell_ref_grid)
  cell_ref_coords = lazy_map(get_node_coordinates,cell_ref_grid)
  cell_ref_conns = lazy_map(get_cell_node_ids,cell_ref_grid)
  cell_type = get_cell_type(ctopo)
  cell_l2g,n_nodes = unstructured_uniform_cell_l2gmap_and_nnodes(ctopo,n,cell_refine_masks)
  n_cells = sum(cell_lncells)
  coords = unstructured_uniform_coordinates(cell_ref_coords,cell_l2g,cell_map,n_nodes)
  conn = unstructured_uniform_connectivity(cell_ref_conns,cell_l2g,n_cells)
  
  f_cell_type = similar(cell_type,n_cells)
  i = 1
  @inbounds for ci in eachindex(cell_lncells)
    lncells = cell_lncells[ci]
    type = cell_type[ci]
    f_cell_type[i:i+lncells-1] .= type
    i += lncells
  end

  c_oriented = OrientationStyle(ctopo) == Oriented()
  ref_oriented = lazy_map(g->OrientationStyle(g)==Oriented(),cell_ref_grid)
  orientation_style = c_oriented && all(ref_oriented) ? Oriented() : NonOriented()

  UnstructuredGridTopology(
    coords,
    conn,
    f_cell_type,
    polytopes,
    orientation_style
  )
end

function unstructured_uniform_refine(cm::DiscreteModel,n::Integer;cell_refine_masks::AbstractVector{Bool})
  @assert length(cell_refine_masks) == num_cells(cm)

  polytopes = get_polytopes(cm)
  # cell_type = get_cell_type(cm)
  ctopo = get_grid_topology(cm)
  cgrid = get_grid(cm)
  cell_map = get_cell_map(cgrid)
  # cell_polytope = CompressedArray(polytopes,cell_type)
  pt_maps = map(cube_simplex_pattern_dimfid,polytopes)
  cmparr_ptrs = collect(get_cell_type(cm))
  cmparr_ptrs[cell_refine_masks] .+= length(polytopes)
  cell_ref_grid = CompressedArray(
    vcat(
      map((p,pt_map)->cube_simplex_reference_grid(p,1,pt_map),polytopes,pt_maps),
      map((p,pt_map)->cube_simplex_reference_grid(p,n,pt_map),polytopes,pt_maps)
    ),
    cmparr_ptrs
  )
  rrules = CompressedArray(
    vcat(
      map((p,g)->RefinementRule(WithoutRefinement(),p,g),polytopes,cell_ref_grid),
      map((p,g)->RefinementRule(GenericRefinement(),p,g),polytopes,cell_ref_grid)
    ),
    cmparr_ptrs
  )

  topo = unstructured_uniform_topology(cell_ref_grid,ctopo,cell_map,n,cell_refine_masks)
  grid = UnstructuredGrid(
    get_vertex_coordinates(topo),
    get_faces(topo,num_cell_dims(topo),0),
    get_reffes(cgrid),
    get_cell_type(topo),
    OrientationStyle(topo)
  )

  glue = blocked_refinement_glue(rrules)
  face_labeling = FaceLabeling(topo)
  model = UnstructuredDiscreteModel(grid,topo,face_labeling)
  AdaptedDiscreteModel(model,cm,glue)
end


# Tester

function test_unstructured_uniform_refinement()
  # _num_*
  n = 3
  polytopes = [SEGMENT,TRI,QUAD,TET,HEX]
  ncells = [3,9,9,27,27]
  nnodes = [4,10,16,20,64]
  n_interior_nodes = [2,1,4,0,8]

  @test all( map(p->_num_cells(p,n),polytopes) .== ncells )
  @test all( map(p->_num_nodes(p,n),polytopes) .== nnodes )
  @test all( map(p->_num_interior_nodes(p,n),polytopes) .== n_interior_nodes )

  @test_throws ErrorException _num_cells(PYRAMID, n)
  @test_throws ErrorException _num_nodes(PYRAMID, n)
  @test_throws ErrorException _num_interior_nodes(PYRAMID, n)

  # points
  xs = VectorValue{2,Float64}[
    (0,0),(0.1,0),(0.1,0.1),(1/3,2/3),(0,1)
  ]
  # 0->-1, 1->0, o->1
  cube_pts = [
    (-1,-1),(1,-1),(1,1),(1,1),(-1,0)
  ]
  simplex_pts = [
    (-1,-1,0),(1,-1,0),(1,1,0),(1,1,1),(-1,0,1)
  ]

  @test xs == sort(xs,lt=_point_isless)
  @test cube_pts == map(_to_cube_pattern,xs)
  @test simplex_pts == map(_to_simplex_pattern,xs)
end