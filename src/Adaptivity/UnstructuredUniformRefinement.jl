# some helper functions

# if `a` is a `CompressedArray`, then replace `b.ptrs` with `a.ptrs`.
_similar(::AbstractArray,b::CompressedArray) = b

function _similar(a::CompressedArray,b::CompressedArray)
  n = Int( length(a.values)/length(b.values) )
  CompressedArray(repeat(b.values,n),a.ptrs)
end

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
############################################################################################
"""
    cube_simplex_pattern_dimfid(p::Polytope{D}) where {D}

Given a simplex or hypercube, return a dictionary to serve as a mapping `f` from 
`pattern` to `(dim, fid)`. For example, for `QUAD`, `f[(-1, -1)]` returns `(0, 1)` because
the origin corresponds to the first 0-dimensional point.
"""
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

@doc raw"""
    cube_simplex_reference_grid(p::Polytope,n::Integer)

Given a simplex or hypercube, return an `n`-partitioned grid, where the vertices
are ordered according to the dimension and id of the face they belong to. 
For example, when n = 2, the result for the `QUAD` and `TRI` are as follows:
3 ---- 6 ---- 4     |     3
|      |      |     |     | \
|      |      |     |     |  \
7 ---- 9 ---- 8     |     5 - 6
|      |      |     |     | \ |\
|      |      |     |     |  \| \
1 ---- 5 ---- 2     |     1 - 4 - 2
      QUAD          |        TRI
"""
function cube_simplex_reference_grid(p::Polytope,n::Integer)
  pt_map = cube_simplex_pattern_dimfid(p)
  atol = 0.25 / n 
  if is_n_cube(p)
    ref_grid = Geometry._compute_linear_grid_from_n_cube(p,n)
    coords = collect(vec(get_node_coordinates(ref_grid)))
    dimfids = map(c->pt_map[_to_cube_pattern(c,atol)],coords)
    conn = Table(get_cell_node_ids(ref_grid))
  else 
    coords,_conn = Geometry._compute_linear_grid_coords_from_simplex(p,n)
    dimfids = map(c->pt_map[_to_simplex_pattern(c,atol)],coords)
    conn = Table(_conn)
  end
  perm = sortperm(broadcast(tuple,dimfids,coords))
  inv_perm = similar(perm)
  inv_perm[perm] = 1:length(perm)
  new_coords = coords[perm]
  new_conn = Table(inv_perm[conn.data],conn.ptrs)
  reffes = [ LagrangianRefFE(Float64,p,1) ]
  orientation_style = NonOriented()
  cell_types = ones(Int8,length(conn))
  has_affine_map = true

  UnstructuredGrid(
    new_coords,
    new_conn,
    reffes,
    cell_types,
    orientation_style;
    has_affine_map
  )
end

"""
    cube_simplex_interior_permutation(p::Polytope,n::Integer)

Given a simplex or hypercube, return the permutation of lids of internal nodes induced by its vertex permutation.
"""
function cube_simplex_interior_permutation(p::Polytope,n::Integer)
  pcoords = get_vertex_coordinates(p)
  reffe = LagrangianRefFE(p)
  shapefuns = get_shapefuns(reffe)
  grid = cube_simplex_reference_grid(p,n)
  coords = get_node_coordinates(grid)[end-_num_interior_nodes(p,n)+1:end]
  map(get_vertex_permutations(p)) do perm
    f = linear_combination(pcoords[perm],shapefuns)
    sortperm(evaluate(f,coords),lt=_point_isless)
  end
end

"""
    compute_d_dface_offsets(ctopo,cell_ref_grid)

Given the coarse topology and cell-wise reference grid, return fine nodal offsets
in each coarse faces.
"""
function compute_d_dface_offsets(ctopo,cell_ref_grid)
  @assert length(cell_ref_grid) == num_cells(ctopo)
  
  Dc = num_cell_dims(ctopo)
  polytopes = get_polytopes(ctopo)
  cell_type = get_cell_type(ctopo) 
  cell_polytope = _similar(cell_ref_grid,CompressedArray(polytopes,cell_type))
  cell_dface_nnodes = lazy_map(cell_ref_grid,cell_polytope) do grid,p
    pt_map = cube_simplex_pattern_dimfid(p)
    coords = get_node_coordinates(grid)
    if is_n_cube(p)
      dimfids = map(c->pt_map[_to_cube_pattern(c)],coords)
    else
      dimfids = map(c->pt_map[_to_simplex_pattern(c)],coords)
    end
    map(d->Int32(count(==((d,1)),dimfids)),0:Dc)::Vector{Int32}
  end
  d_dface_offsets = Vector{Vector{Int32}}(undef,Dc)
  for d in 1:Dc
    (;data,ptrs) = Table(get_faces(ctopo,d,Dc))
    d_offsets = similar(ptrs,Int32)
    d_offsets[1] = d == 1 ? num_vertices(ctopo)+1 : d_dface_offsets[d-1][end]
    @inbounds for dfi = 1:length(d_offsets)-1
      I = findfirst(view(data,ptrs[dfi]:ptrs[dfi+1]-1)) do i 
        num_vertices(cell_polytope[i]) < num_nodes(cell_ref_grid[i])
      end
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

"""
    compute_cell_offsets(ctopo,cell_ref_grid)

Given the coarse topology and cell-wise reference grid, return cell-wise fine nodal offsets.
"""
function compute_cell_offsets(ctopo,cell_ref_grid)
  @assert length(cell_ref_grid) == num_cells(ctopo)

  Dc = num_cell_dims(ctopo)
  cell_nnodes = lazy_map(num_nodes,cell_ref_grid)
  offsets = Vector{Int32}(undef,num_faces(ctopo,Dc)+1)
  offsets[1] = 1
  @inbounds for ci in 1:length(offsets)-1
    offsets[ci+1] = offsets[ci] + cell_nnodes[ci]
  end
  offsets
end

"""
    unstructured_refine_cell_l2gmap_and_nnodes(
      ctopo,
      cell_ref_grid,
      cell_dface_permutations::AbstractVector)

Cell-wise local-to-global fine nodal mappings and the number of all fine nodes.
"""
function unstructured_refine_cell_l2gmap_and_nnodes(
  ctopo,
  cell_ref_grid,
  cell_dface_permutations::AbstractVector)
  @assert num_cells(ctopo) == length(cell_dface_permutations) == length(cell_ref_grid)

  Dc = num_cell_dims(ctopo)
  d_df_goffsets = compute_d_dface_offsets(ctopo,cell_ref_grid)
  l2g_ptrs = compute_cell_offsets(ctopo,cell_ref_grid)
  l2g_data = Vector{Int32}(undef,l2g_ptrs[end]-1)
  d_c2df = ntuple(d->Table(get_faces(ctopo,Dc,d)),Val{Dc}())
  d_c2perm = ntuple(d->Table(get_cell_permutations(ctopo,d)),Val{Dc}())
  c2n = Table(get_faces(ctopo,Dc,0))
  cell_polytope = CompressedArray(get_polytopes(ctopo),get_cell_type(ctopo))

  @inbounds Threads.@threads for ci in 1:length(l2g_ptrs)-1
    lo = l2g_ptrs[ci]
    # vertices
    pini = c2n.ptrs[ci]
    l = c2n.ptrs[ci+1]-pini
    copyto!(l2g_data,lo,c2n.data,pini,l)
    lo += l

    # 1-face ~ (Dc)-face
    dface_perm = cell_dface_permutations[ci]
    if num_vertices(cell_polytope[ci]) < num_nodes(cell_ref_grid[ci])
      for d in 1:Dc
        c2perm = d_c2perm[d]
        c2df = d_c2df[d]
        goffsets = d_df_goffsets[d]
        pini = c2df.ptrs[ci]
        ndf = c2df.ptrs[ci+1] - pini
        for i in 0:ndf-1
          fi = c2df.data[pini+i]
          tp = c2perm.data[pini+i]
          perm = dface_perm[d][tp]
          gini = goffsets[fi]
          nn = goffsets[fi+1] - gini
          for j in 0:nn-1
            l2g_data[lo+perm[j+1]-1] = gini + j
          end
          lo += nn
        end
      end
    end
  end

  n_nodes = d_df_goffsets[end][end] - 1
  Table(l2g_data,l2g_ptrs),n_nodes
end

function unstructured_refine_coordinates(cell_ref_coords,cell_l2g,cell_map,n_nodes)
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

function unstructured_refine_connectivity(cell_ref_conns,cell_l2g,n_cells)
  cell_ref_conns = lazy_map(Table, cell_ref_conns)
  conn_ptrs = Vector{Int32}(undef,n_cells+1)
  conn_ptrs[1] = 1
  fci = 1
  @inbounds for ci in eachindex(cell_ref_conns)
    (;ptrs) = cell_ref_conns[ci]
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

function unstructured_refine_topology(
  cell_ref_grid,
  ctopo,
  cell_map,
  cell_dface_permutations::AbstractVector)
  polytopes = get_polytopes(ctopo)
  @notimplementedif !all( map(p->is_n_cube(p) || is_simplex(p),polytopes) )

  cell_lncells = lazy_map(num_cells,cell_ref_grid)
  cell_ref_coords = lazy_map(get_node_coordinates,cell_ref_grid)
  cell_ref_conns = lazy_map(get_cell_node_ids,cell_ref_grid)
  cell_type = get_cell_type(ctopo)
  cell_l2g,n_nodes = unstructured_refine_cell_l2gmap_and_nnodes(
    ctopo,
    cell_ref_grid,
    cell_dface_permutations
  )
  n_cells = sum(cell_lncells)
  coords = unstructured_refine_coordinates(cell_ref_coords,cell_l2g,cell_map,n_nodes)
  conn = unstructured_refine_connectivity(cell_ref_conns,cell_l2g,n_cells)
  
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

"""
    unstructured_refine(
      cm::DiscreteModel,
      cell_ref_grid::AbstractVector{<:UnstructuredGrid},
      cell_dface_permutations::AbstractVector;
      has_affine_map::Union{Nothing,Bool}=nothing)

Normally, given the cell-wise reference grids, a refined mesh can be constructed via `get_cell_map(cm)`. 
However, since the `get_cell_map(cm)` does not contain permutation information for faces, we need 
`cell_dface_permutations` to remap the newly generated points within each face accordingly. This approache
may lead to the presence of hanging nodes. However, this property is preserved as it might be useful 
in future applications. 

Note: if a face belongs to cells marked for refinement, it is always assumed that the newly introduced points
within the face also belong to these refined cells. 
"""
function unstructured_refine(
  cm::DiscreteModel,
  cell_ref_grid::AbstractVector{<:UnstructuredGrid},
  cell_dface_permutations::AbstractVector;
  has_affine_map::Union{Nothing,Bool}=nothing)
  @assert num_cells(cm) == length(cell_dface_permutations) == length(cell_ref_grid)

  ctopo = get_grid_topology(cm)
  cgrid = get_grid(cm)
  cell_map = get_cell_map(cgrid)
  _is_affine(fs) = isconcretetype(typeof(fs)) && fs isa AbstractArray{<:AffineField}
  cell_polytope = _similar(
    cell_ref_grid,
    CompressedArray(get_polytopes(ctopo),get_cell_type(ctopo))
  )

  rrules = lazy_map(cell_polytope,cell_ref_grid) do p,g
    if num_vertices(p) < num_nodes(g)
      RefinementRule(GenericRefinement(),p,g)
    else
      RefinementRule(WithoutRefinement(),p,g)
    end
  end
  
  topo = unstructured_refine_topology(
    cell_ref_grid,
    ctopo,
    cell_map,
    cell_dface_permutations
  )
  if has_affine_map isa Nothing
    cell_ref_is_affine = lazy_map(g->_is_affine(get_cell_map(g)),cell_ref_grid)
    has_affine_map = all(cell_ref_is_affine) && _is_affine(cell_map)
  end

  grid = UnstructuredGrid(
    get_vertex_coordinates(topo),
    get_faces(topo,num_cell_dims(topo),0),
    get_reffes(cgrid),
    get_cell_type(topo),
    OrientationStyle(topo);
    has_affine_map
  )

  glue = blocked_refinement_glue(rrules)
  face_labeling = FaceLabeling(topo)
  model = UnstructuredDiscreteModel(grid,topo,face_labeling)
  AdaptedDiscreteModel(model,cm,glue)
end


"""
    unstructured_uniform_refine(
      cm::DiscreteModel,
      n::Integer,
      cell_refine_masks::AbstractVector{Bool};
      kwargs...)

Uniformly refine the given discrete model into `n` parts per dimension. It determine which coarse cells 
need to be refined based on `cell_refine_masks`.
"""
function unstructured_uniform_refine(
  cm::DiscreteModel,
  n::Integer,
  cell_refine_masks::AbstractVector{Bool};
  kwargs...)

  polytopes = get_polytopes(cm)
  cell_type = get_cell_type(cm)
  cmparr_ptrs = collect(cell_type)
  cmparr_ptrs[cell_refine_masks] .+= length(polytopes)
  cell_ref_grid = CompressedArray(
    vcat(
      map(p->cube_simplex_reference_grid(p,1),polytopes),
      map(p->cube_simplex_reference_grid(p,n),polytopes)
    ),
    cmparr_ptrs
  )
  dface_permutations = map(polytopes) do p 
    map(fp->(cube_simplex_interior_permutation(fp,n)),get_reffaces(p)[2:end])
  end
  cell_dface_permutations = CompressedArray(dface_permutations,cell_type)
  unstructured_refine(cm,cell_ref_grid,cell_dface_permutations;kwargs...)
end

# Tester

function test_unstructured_uniform_refinement()
  # _similar
  A = CompressedArray([10,11,12,13],[1,2,3,4])
  B = CompressedArray([20,21],[1,2])
  C = [30,31,32]
  @test _similar(C,A) == A
  @test _similar(A,B) == [20,21,20,21]
  @test_throws InexactError _similar(B,A)

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