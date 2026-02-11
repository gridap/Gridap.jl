# High level interfaces

function refine(model::CartesianDiscreteModel{Dc}, cell_partition::Int=2) where Dc
  partition = Tuple(fill(cell_partition,Dc))
  return refine(model,partition)
end

function refine(model::CartesianDiscreteModel{Dc}, cell_partition::Tuple) where Dc
  desc = Geometry.get_cartesian_descriptor(model)
  nC   = desc.partition

  # Refinement Glue
  f2c_cell_map, fcell_to_child_id = _create_cartesian_f2c_maps(nC,cell_partition)
  faces_map = [(d==Dc) ? f2c_cell_map : Int[] for d in 0:Dc]
  reffe     = LagrangianRefFE(Float64,first(get_polytopes(model)),1)
  rrules    = RefinementRule(reffe,cell_partition)
  glue = AdaptivityGlue(faces_map,fcell_to_child_id,rrules)

  # Refined model
  domain     = _get_cartesian_domain(desc)
  _model_ref = CartesianDiscreteModel(domain,cell_partition.*nC)

  # Propagate face labels
  coarse_labels = get_face_labeling(model)
  coarse_topo   = get_grid_topology(model)
  fine_topo     = get_grid_topology(_model_ref)
  fine_labels   = refine_face_labeling(coarse_labels,glue,coarse_topo,fine_topo)

  model_ref = CartesianDiscreteModel(get_grid(_model_ref),fine_topo,fine_labels)
  return AdaptedDiscreteModel(model_ref,model,glue)
end

function refine(model::UnstructuredDiscreteModel,cell_partition::Int)
  @check cell_partition >= 1 
  if cell_partition == 1
    model
  else
    cell_refine_masks = Fill(true,num_cells(model))
    uniformly_refine(model,cell_partition,cell_refine_masks)
  end
end

function refine(model::Geometry.DiscreteModelMock,cell_partition::Int)
  @check cell_partition >= 1
  if cell_partition == 1
    model
  else
    cell_refine_masks = Fill(true,num_cells(model))
    uniformly_refine(model,cell_partition,cell_refine_masks)
  end
end

"""
  uniformly_refine(
    cm::DiscreteModel,
    n::Integer,
    cell_refine_masks::AbstractVector{T};
    has_affine_map::Union{Nothing, Bool} = nothing
  ) where T <: Union{Bool, Int}

Uniformly refine the cells of the discrete model `cm` marked by `cell_refine_masks` `n` times.
If `has_affine_map` is not provided, it is automatically determined.
"""
function uniformly_refine(
  cm::DiscreteModel,
  n::Integer,
  cell_refine_masks::AbstractVector{T};
  has_affine_map::Union{Nothing, Bool} = nothing) where T <: Union{Bool, Int}

  @check n >= 1 "The number of uniform refinements must be at least one."
  if n == 1
    return cm
  end
  ctopo = get_grid_topology(cm)
  cgrid = get_grid(cm)
  cell_map = get_cell_map(cgrid)
  polytopes = get_polytopes(ctopo)
  cmparr_ptrs = collect(get_cell_type(cm))
  cmparr_ptrs[cell_refine_masks] .+= length(polytopes)

  without_rr = map(polytopes) do p
    RefinementRule(WithoutRefinement(), p, compute_reference_grid(p, 1))
  end
  generic_rr = map(polytopes) do p
    RefinementRule(GenericRefinement(), p, compute_reference_grid(p, n))
  end
  cell_rr = CompressedArray(vcat(without_rr, generic_rr), cmparr_ptrs)
  _is_affine(fs) = isconcretetype(typeof(fs)) && fs isa AbstractArray{<:AffineField}
  glue = blocked_refinement_glue(cell_rr)

  if has_affine_map isa Nothing
    cell_ref_is_affine = lazy_map(rr->_is_affine(get_cell_map(rr)), cell_rr)
    has_affine_map = all(cell_ref_is_affine) && _is_affine(cell_map)
  end

  grid, topo = _uniformly_refine_grid_topology(
    cell_rr,
    ctopo,
    cgrid,
    cell_map,
    has_affine_map,
  )
  clabeling = get_face_labeling(cm)
  labeling = refine_face_labeling(clabeling, glue, ctopo, topo)
  model = UnstructuredDiscreteModel(grid, topo, labeling)
  AdaptedDiscreteModel(model, cm, glue)
end

# Implementation details

function _get_cartesian_domain(desc::CartesianDescriptor{D}) where D
  origin = desc.origin
  corner = origin + VectorValue(desc.sizes .* desc.partition)
  domain = Vector{eltype(origin)}(undef,2*D)
  for d in 1:D
    domain[d*2-1] = origin[d]
    domain[d*2]   = corner[d]
  end
  return Tuple(domain)
end

@generated function _c2v(idx::Union{NTuple{N,T},CartesianIndex{N}},sizes::NTuple{N,T}) where {N,T}
  res = :(idx[1])
  for d in 1:N-1
    ik = :((idx[$(d+1)]-1))
    for k in 1:d
        ik = :($ik * sizes[$k])
    end
    res = :($res + $ik)
  end
  return res
end

@generated function _create_cartesian_f2c_maps(nC::NTuple{N,T},ref::NTuple{N,T}) where {N,T}
  J_f2c   = Meta.parse(prod(["(",["1+(I[$k]-1)÷ref[$k]," for k in 1:N]...,")"]))
  J_child = Meta.parse(prod(["(",["1+(I[$k]-1)%ref[$k]," for k in 1:N]...,")"]))

  return :(begin
    nF = nC .* ref
    f2c_map   = Vector{Int}(undef,prod(nF))
    child_map = Vector{Int}(undef,prod(nF))

    for (i,I) in enumerate(CartesianIndices(nF))
      J_f2c   = $J_f2c
      J_child = $J_child
      f2c_map[i] = _c2v(J_f2c,nC)
      child_map[i] = _c2v(J_child,ref)
    end

    return f2c_map, child_map
  end)
end

# Return the local indices of the fine points within each face.
function _uniform_d_dface_own_lid(rr)
  # NOTE: The nodes on the faces in each dimension 
  # need to be ordered in a consistent manner.
  grid = get_ref_grid(rr)
  Dc = num_cell_dims(grid)
  atol = 0.1/(2*num_nodes(grid))^(1/Dc)
  function _point_isless(x, y)
    @inbounds for i ∈ Dc:-1:1
      if x[i] < y[i]-atol
        return true
      elseif x[i] > y[i]+atol
        return false
      end
    end
    true
  end

  face_vertices = get_face_vertices(rr)
  face_coords = get_face_coordinates(rr)
  face_lids = Vector{Vector{Int32}}(undef, length(face_vertices))
  seen = Int32[]
  @inbounds for (i, fvs) in enumerate(face_vertices)
    lids = Int32[]
    perm = sortperm(face_coords[i], lt = _point_isless)
    for (j, v) in enumerate(fvs)
      if !(v ∈ seen)
        push!(lids, perm[j])
        push!(seen, v)
      end
    end
    face_lids[i] = lids
  end
  face_lids
end

# For each dimension and each d-face, return the offsets
# of the global indices associated with the refined mesh.
function _uniform_d_dface_new_offsets(
  ctopo,
  cell_rr,
  cell_dface_n_innodes)

  @assert length(cell_rr) == num_cells(ctopo)
  Dc = num_cell_dims(ctopo)
  d_dface_offsets = Vector{Vector{Int32}}(undef, Dc)
  for d in 1:Dc
    (; data, ptrs) = Table(get_faces(ctopo, d, Dc))
    d_offsets = similar(ptrs, Int32)
    d_offsets[1] = d==1 ? num_vertices(ctopo)+1 : d_dface_offsets[d-1][end]
    @inbounds for dfi ∈ 1:(length(d_offsets)-1)
      I = findfirst(view(data, ptrs[dfi]:(ptrs[dfi+1]-1))) do i
        cell_rr[i].T isa GenericRefinement
      end
      if !isnothing(I)
        ci = data[ptrs[dfi]+I-1]
        d_offsets[dfi+1] = d_offsets[dfi] + cell_dface_n_innodes[ci][d+1]
      else
        d_offsets[dfi+1] = d_offsets[dfi]
      end
    end
    d_dface_offsets[d] = d_offsets
  end
  d_dface_offsets
end

# Return the `ptrs` to construct the Table `cell_l2g_map`.
function _uniform_l2g_ptrs(ctopo, cell_rr)
  @assert length(cell_rr) == num_cells(ctopo)
  cell_nnodes = lazy_map(r->num_nodes(get_ref_grid(r)), cell_rr)
  offsets = Vector{Int32}(undef, num_cells(ctopo)+1)
  offsets[1] = 1
  @inbounds for ci in 1:(length(offsets)-1)
    offsets[ci+1] = offsets[ci] + cell_nnodes[ci]
  end
  offsets
end

# Return the local-to-global map and the number of nodes
# of the uniformly refined mesh.
function _uniform_cell_l2gmap_and_nnodes(
  ctopo,
  cell_rr,
  cell_face_own_vertex_permutations)

  @assert num_cells(ctopo) == length(cell_rr)
  Dc = num_cell_dims(ctopo)
  cell_dface_n_innodes = lazy_map(cell_face_own_vertex_permutations, cell_rr) do fperms, rr
    poly = get_polytope(rr)
    offsets = get_offsets(poly) .+ 1
    map(o -> length(first(fperms[o])), offsets)
  end
  cell_dim_ranges = lazy_map(r->get_dimranges(get_polytope(r)), cell_rr)
  d_df_goffsets = _uniform_d_dface_new_offsets(ctopo, cell_rr, cell_dface_n_innodes)
  l2g_ptrs = _uniform_l2g_ptrs(ctopo, cell_rr)
  l2g_data = Vector{Int32}(undef, l2g_ptrs[end]-1)
  d_c2df = ntuple(d->Table(get_faces(ctopo, Dc, d)), Val{Dc}())
  d_c2perm = ntuple(d->Table(get_cell_permutations(ctopo, d)), Val{Dc}())
  c2n = Table(get_faces(ctopo, Dc, 0))

  @inbounds Threads.@threads for ci in 1:(length(l2g_ptrs)-1)
    face_own_vertex_perm = cell_face_own_vertex_permutations[ci]
    # vertices
    for i in cell_dim_ranges[ci][1]
      lo = face_own_vertex_perm[i][1][1]
      l2g_data[l2g_ptrs[ci]+lo-1] = c2n.data[c2n.ptrs[ci]+i-1]
    end

    # 1-Dc faces
    for d ∈ 1:Dc
      rng = cell_dim_ranges[ci][d+1]
      df_perm = d_c2perm[d]
      pini = df_perm.ptrs[ci]
      c2df = d_c2df[d]
      for (lfi, cfi) in enumerate(rng)
        gfi = c2df.data[c2df.ptrs[ci]+lfi-1]
        go = d_df_goffsets[d][gfi]
        p = df_perm.data[pini+lfi-1]
        for (lo, lni) in enumerate(face_own_vertex_perm[cfi][p])
          l2g_data[l2g_ptrs[ci]+lni-1] = go + lo - 1
        end
      end
    end
  end

  n_nodes = d_df_goffsets[end][end] - 1
  Table(l2g_data, l2g_ptrs), n_nodes
end

# Return the coordinates of the uniformly refined mesh.
function _uniform_coordinates(cell_ref_coords, cell_l2g, cell_map, n_nodes)
  coords = similar(first(cell_ref_coords), n_nodes)
  cell_l2g = Table(cell_l2g)
  @inbounds Threads.@threads for ci in eachindex(cell_l2g)
    ref_co = cell_ref_coords[ci]
    (; data, ptrs) = cell_l2g
    f = cell_map[ci]
    pini = ptrs[ci]
    l = ptrs[ci+1]-pini
    for j in 1:l
      coords[data[pini+j-1]] = f(ref_co[j])
    end
  end
  coords
end

# Return the connectivity of the uniformly refined mesh.
function _uniform_connectivity(cell_ref_conns, cell_l2g, n_cells)
  cell_ref_conns = lazy_map(Table, cell_ref_conns)
  conn_ptrs = Vector{Int32}(undef, n_cells+1)
  conn_ptrs[1] = 1
  fci = 1
  @inbounds for ci in eachindex(cell_ref_conns)
    (; ptrs) = cell_ref_conns[ci]
    for i in 1:(length(ptrs)-1)
      conn_ptrs[fci+1] = conn_ptrs[fci] + ptrs[i+1] - ptrs[i]
      fci += 1
    end
  end

  conn_data = Vector{Int32}(undef, conn_ptrs[end]-1)
  cache = array_cache(cell_l2g)
  go = 1
  @inbounds for ci in eachindex(cell_ref_conns)
    (; data) = cell_ref_conns[ci]
    n = length(data)
    l2g = getindex!(cache, cell_l2g, ci)
    conn_data[go:(go+n-1)] = l2g[data]
    go += n
  end
  Table(conn_data, conn_ptrs)
end

# Construct the grid and topology for the uniformly refined grid.
function _uniformly_refine_grid_topology(
  cell_rr,
  ctopo,
  cgrid,
  cell_map,
  has_affine_map)

  cell_ref_grid = lazy_map(get_ref_grid, cell_rr)
  cell_lncells = lazy_map(num_cells, cell_ref_grid)
  cell_ref_coords = lazy_map(get_node_coordinates, cell_ref_grid)
  cell_ref_conns = lazy_map(get_cell_node_ids, cell_ref_grid)
  cell_type = get_cell_type(ctopo)
  cell_face_own_vertex_permutations = lazy_map(cell_rr) do rr
    face_vertices = get_face_vertices(rr)
    lids = _uniform_d_dface_own_lid(rr)
    vertex_perm = get_face_vertex_permutations(rr)
    own_vertex_perm = similar(vertex_perm)
    @inbounds for i ∈ eachindex(vertex_perm)
      perms = vertex_perm[i]
      lid = lids[i]
      vertices = face_vertices[i]
      own_vertex_perm[i] = map(p -> vertices[p[lid]], perms)
    end
    own_vertex_perm
  end

  cell_l2g, n_nodes = _uniform_cell_l2gmap_and_nnodes(
    ctopo,
    cell_rr,
    cell_face_own_vertex_permutations,
  )
  n_cells = sum(cell_lncells)
  coords = _uniform_coordinates(cell_ref_coords, cell_l2g, cell_map, n_nodes)
  conn = _uniform_connectivity(cell_ref_conns, cell_l2g, n_cells)
  f_cell_type = similar(cell_type, n_cells)
  i = 1
  @inbounds for ci in eachindex(cell_lncells)
    lncells = cell_lncells[ci]
    type = cell_type[ci]
    f_cell_type[i:(i+lncells-1)] .= type
    i += lncells
  end
  c_oriented = OrientationStyle(ctopo) == Oriented()
  ref_oriented = lazy_map(g->OrientationStyle(g)==Oriented(), cell_ref_grid)
  orientation_style = c_oriented && all(ref_oriented) ? Oriented() : NonOriented()

  ftopo = UnstructuredGridTopology(
    coords,
    conn,
    f_cell_type,
    get_polytopes(ctopo),
    orientation_style,
  )

  # handle periodic cases
  if num_vertices(ctopo) != num_nodes(cgrid)
    grid_cell_l2g, grid_n_nodes = _uniform_cell_l2gmap_and_nnodes(
      UnstructuredGridTopology(cgrid),
      cell_rr,
      cell_face_own_vertex_permutations,
    )
    grid_coords = _uniform_coordinates(cell_ref_coords, grid_cell_l2g, cell_map, grid_n_nodes)
    grid_conn = _uniform_connectivity(cell_ref_conns, grid_cell_l2g, n_cells)
  else
    grid_coords = get_vertex_coordinates(ftopo)
    grid_conn = get_faces(ftopo, num_cell_dims(ftopo), 0)
  end
  fgrid = UnstructuredGrid(
    grid_coords,
    grid_conn,
    get_reffes(cgrid),
    f_cell_type,
    orientation_style;
    has_affine_map,
  )

  return fgrid, ftopo
end