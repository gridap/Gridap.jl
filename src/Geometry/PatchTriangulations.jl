
"""
    struct PatchTopology{Dc,Dp} <: GridapType
      topo :: GridTopology{Dc,Dp}
      d_to_patch_to_dfaces :: Vector{Table{Int32,Vector{Int32},Vector{Int32}}}
    end

# Fields: 

- `topo`: Underlying grid topology
- `d_to_patch_to_dfaces`: For each dimension `d`, a table mapping each patch to it's d-faces.
"""
struct PatchTopology{Dc,Dp,A} <: GridapType
  topo :: GridTopology{Dc,Dp}
  d_to_patch_to_dfaces :: Vector{Table{Int32,Vector{Int32},Vector{Int32}}}
  metadata :: A
  function PatchTopology(
    topo::GridTopology{Dc,Dp},
    d_to_patch_to_dfaces::Vector{Table{Int32,Vector{Int32},Vector{Int32}}},
    metadata=nothing
  ) where {Dc,Dp}
    A = typeof(metadata)
    new{Dc,Dp,A}(topo,d_to_patch_to_dfaces,metadata)
  end
end

function PatchTopology(topo::GridTopology{Dc},patch_cells::Table,metadata=nothing) where Dc
  d_to_patch_to_dpfaces = Vector{Table{Int32,Vector{Int32},Vector{Int32}}}(undef,Dc+1)
  d_to_patch_to_dpfaces[Dc+1] = patch_cells
  PatchTopology(topo,d_to_patch_to_dpfaces,metadata)
end

struct StarPatchMetadata
  root_dim :: Int8
  patch_roots :: Vector{Int32}
end

function PatchTopology(
  ::Type{ReferenceFE{Df}},model::DiscreteModel;
  labels = get_face_labeling(model), tags = nothing
) where Df
  Dc = num_cell_dims(model)
  topo = get_grid_topology(model)
  patch_cells = get_faces(topo,Df,Dc)
  patch_roots = collect(Int32,1:num_faces(topo,Df))
  if !isnothing(tags)
    patches = findall(get_face_mask(labels,tags,Df))
    patch_cells = patch_cells[patches]
    patch_roots = patch_roots[patches]
  end
  metadata = StarPatchMetadata(Int8(Df),patch_roots)
  PatchTopology(topo,patch_cells,metadata)
end

function PatchTopology(model::DiscreteModel;kwargs...)
  D = num_cell_dims(model)
  PatchTopology(ReferenceFE{D},model;kwargs...)
end

function InterfacePatchTopology(model::DiscreteModel; kwargs...)
  ptopo = PatchTopology(model; kwargs...)
  InterfacePatchTopology(ptopo)
end

function InterfacePatchTopology(ptopo::PatchTopology; patch_ids = 1:num_patches(ptopo))
  topo = ptopo.topo
  D = num_cell_dims(ptopo)
  face_to_cells = get_faces(topo,D-1,D)

  patch_faces = get_patch_faces(ptopo,D-1)
  pface_to_isboundary, pface_to_lcell = get_patch_boundary_info(ptopo)
  face_to_patches = Arrays.inverse_table(patch_faces,num_faces(topo,D-1))
  pface_to_face = patch_faces.data

  ni = 0
  id_to_iface = Dict{UInt,Int32}()
  pface_to_iface = zeros(Int32, num_faces(ptopo,D-1))
  pface_to_cell = zeros(Int32, num_faces(ptopo,D-1))
  for patch in patch_ids
    for pface in datarange(patch_faces,patch)
      if !pface_to_isboundary[pface]
        continue
      end
      face = pface_to_face[pface]
      patches = view(face_to_patches,face)
      cells = view(face_to_cells,face)
      key = hash(Set(patches))
      iface = get!(id_to_iface, key, ni+1)
      pface_to_iface[pface] = iface
      pface_to_cell[pface] = cells[pface_to_lcell[pface]]
      ni += isequal(iface,ni+1)
    end
    empty!(id_to_iface)
  end
  iface_to_pfaces = Arrays.inverse_table(pface_to_iface,ni)

  d_to_iface_to_dfaces = Vector{Table{Int32,Vector{Int32},Vector{Int32}}}(undef,D+1)
  d_to_iface_to_dfaces[D+1] = Table(lazy_map(unique,lazy_map(Broadcasting(Reindex(pface_to_cell)), iface_to_pfaces)))
  d_to_iface_to_dfaces[D] = Table(lazy_map(Broadcasting(Reindex(pface_to_face)), iface_to_pfaces))
  
  return PatchTopology(topo,d_to_iface_to_dfaces,nothing)
end

function num_patches(ptopo::PatchTopology)
  length(get_patch_cells(ptopo))
end

num_cell_dims(::PatchTopology{Dc}) where Dc = Dc
num_point_dims(::PatchTopology{Dc,Dp}) where {Dc,Dp} = Dp

num_cells(ptopo::PatchTopology{Dc}) where Dc = num_faces(ptopo,Dc)
num_faces(ptopo::PatchTopology,d) = length(get_patch_faces(ptopo,d).data)

get_patch_cells(ptopo::PatchTopology{Dc}) where Dc = get_patch_faces(ptopo,Dc)
get_patch_facets(ptopo::PatchTopology{Dc}) where Dc = get_patch_faces(ptopo,Dc-1)

"""
    get_patch_faces(ptopo::PatchTopology,d::Integer) -> patch_faces
"""
function get_patch_faces(ptopo::PatchTopology,d::Integer)
  if isassigned(ptopo.d_to_patch_to_dfaces,d+1)
    patch_faces = ptopo.d_to_patch_to_dfaces[d+1]
  else
    patch_faces = generate_patch_faces(ptopo,d)
    ptopo.d_to_patch_to_dfaces[d+1] = patch_faces
  end
  return patch_faces
end

"""
    patch_reindex(ptopo::PatchTopology{Dc},face_to_data,Df=Dc) -> pface_to_data
"""
function patch_reindex(ptopo::PatchTopology{Dc},face_to_data,Df=Dc) where Dc
  patch_faces = get_patch_faces(ptopo,Df)
  pface_to_data = lazy_map(Reindex(face_to_data),patch_faces.data)
  return pface_to_data
end

"""
    patch_extend(PD::PatchTopology{Dc},patch_to_data,Df=Dc) -> pface_to_data
"""
function patch_extend(ptopo::PatchTopology{Dc},patch_to_data,Df=Dc) where Dc
  pface_to_patch = Arrays.block_identity_array(get_patch_faces(ptopo,Df).ptrs)
  pface_to_data = lazy_map(Reindex(patch_to_data),pface_to_patch)
  return pface_to_data
end

"""
    is_cover(ptopo::PatchTopology)
    is_cover(topo::GridTopology, patch_cells::Table)

Returns `true` if the given patch topology is a cover of the underlying topology, i.e 
if every cell in the topology is contained in at least one patch.
"""
function is_cover(ptopo::PatchTopology)
  patch_cells = Geometry.get_patch_cells(ptopo)
  is_cover(ptopo.topo, patch_cells)
end

function is_cover(topo::GridTopology, patch_cells)
  cache = array_cache(patch_cells)
  cell_is_covered = fill(false, num_cells(topo))
  for patch in eachindex(patch_cells)
    cells = getindex!(cache, patch_cells, patch)
    cell_is_covered[cells] .= true
  end
  return all(cell_is_covered)
end

"""
    is_disjoint(ptopo::PatchTopology)
    is_disjoint(topo::GridTopology, patch_cells::Table)

Check if the given patch topology is disjoint, i.e if no cell in the topology
is contained in more than one patch.
"""
function is_disjoint(ptopo::PatchTopology)
  patch_cells = Geometry.get_patch_cells(ptopo)
  is_disjoint(ptopo.topo, patch_cells)
end

function is_disjoint(topo::GridTopology, patch_cells)
  cache = array_cache(patch_cells)
  cell_is_covered = fill(false, num_cells(topo))
  for patch in eachindex(patch_cells)
    cells = getindex!(cache, patch_cells, patch)
    if any(view(cell_is_covered,cells))
      return false
    end
    cell_is_covered[cells] .= true
  end
  return true
end

"""
  is_partition(ptopo::PatchTopology; kwargs...)
  is_partition(topo::GridTopology, patch_cells::Table; fail_fast = true)

Check if the given patch topology is a valid partition of the underlying topology.

To be a valid partition, the patches

  - must cover the whole topology
  - must be disjoint (i.e non-overlapping)
  - must be connected (a patch cannot be split)

If `fail_fast` is `true`, the function will exit as soon as it finds a patch that is not connected.
Otherwise, it will check all patches and print a warning with the indices of the bad patches.

"""
function is_partition(ptopo::PatchTopology; kwargs...)
  patch_cells = Geometry.get_patch_cells(ptopo)
  is_partition(ptopo.topo, patch_cells; kwargs...)
end

function is_partition(topo::GridTopology, patch_cells; fail_fast = true, require_cover = true)
  if (require_cover && !is_cover(topo, patch_cells)) || !is_disjoint(topo, patch_cells)
    return false # Not disjoint
  end

  # Check patch connectivity
  # We could return false at the first bad patch, but I'd rather get 
  # all bad patch indices for debugging purposes.
  D = num_cell_dims(topo)
  cell_to_facets = get_faces(topo,D,D-1)
  facet_to_cells = get_faces(topo,D-1,D)

  cache = array_cache(patch_cells)
  bad_patches = Int[]
  for patch in eachindex(patch_cells)
    cells = getindex!(cache, patch_cells, patch)
    G = Graph(compute_graph(Set(cells), cell_to_facets, facet_to_cells))
    if !is_connected(G)
      push!(bad_patches, patch)
      fail_fast && return false # Stop at the first bad patch
    end
  end

  if !isempty(bad_patches)
    @warn "The following patches are not connected: $(bad_patches)"
    return false
  end
  return true
end

function generate_patch_faces(ptopo::PatchTopology{Dc},Df) where Dc
  cell_to_faces = get_faces(ptopo.topo,Dc,Df)
  patch_cells = get_patch_cells(ptopo)
  return Arrays.merge_entries(cell_to_faces, patch_cells; acc = SortedSet{Int32}())
end

function generate_patch_faces(ptopo::PatchTopology,dimfrom::Integer,dimto::Integer)
  npatches = num_patches(ptopo)

  patch_to_fface = get_patch_faces(ptopo,dimfrom)
  patch_to_tface = get_patch_faces(ptopo,dimto)
  fface_to_tface = get_faces(ptopo.topo,dimfrom,dimto)

  k = 1
  ptrs = zeros(Int32, num_faces(ptopo,dimfrom)+1)
  for patch in 1:npatches
    ffaces = view(patch_to_fface,patch)
    tfaces = view(patch_to_tface,patch)
    for fface in ffaces
      f_tfaces = view(fface_to_tface,fface)
      for f in f_tfaces
        pos = searchsortedfirst(tfaces,f)
        ptrs[k+1] += Int32(pos != length(tfaces)+1)
      end
      k += 1
    end
  end
  @assert k == length(ptrs)
  Arrays.length_to_ptrs!(ptrs)

  k = 1
  data = Vector{Int32}(undef, ptrs[end])
  for patch in 1:npatches
    ffaces = view(patch_to_fface,patch)
    tfaces = view(patch_to_tface,patch)
    offset = patch_to_tface.ptrs[patch] - 1
    for fface in ffaces
      f_tfaces = view(fface_to_tface,fface)
      for f in f_tfaces
        pos = searchsortedfirst(tfaces,f)
        if pos != length(tfaces)+1
          data[ptrs[k]] = pos + offset
          ptrs[k] += 1
        end
      end
      k += 1
    end
  end
  @assert k == length(ptrs)
  Arrays.rewind_ptrs!(ptrs)

  return Table(data,ptrs)
end

function compute_isboundary_face(ptopo::PatchTopology{Dc},d::Integer) where Dc
  # Cells
  if d == Dc 
    return fill(false,num_faces(ptopo,Dc))
  end

  # Facets
  isboundary_facet, _ = get_patch_boundary_info(ptopo)
  if d == Dc-1 
    return isboundary_facet
  end
  
  # Faces: Boundary if belongs to a boundary facet
  facet_to_faces = generate_patch_faces(ptopo,Dc-1,d)
  face_mask = fill(false,num_faces(ptopo,d))
  for (facet,mask) in enumerate(isboundary_facet)
    if mask
      faces = view(facet_to_faces,facet)
      face_mask[faces] .= true
    end
  end
  return face_mask
end

function get_patch_boundary_info(ptopo::PatchTopology{Dc}) where Dc
  patch_cells = get_patch_cells(ptopo)
  patch_faces = get_patch_faces(ptopo,Dc-1)
  face_to_cells = get_faces(ptopo.topo,Dc-1,Dc)

  n_pfaces = num_faces(ptopo,Dc-1)
  pface_to_isboundary = fill(false,n_pfaces)
  pface_to_lcell = fill(Int8(1),n_pfaces)

  k = 1
  for patch in 1:num_patches(ptopo)
    pcells = view(patch_cells,patch)
    pfaces = view(patch_faces,patch)
    for pface in pfaces
      pface_cells = view(face_to_cells,pface)
      lcells = findall(c -> c in pcells, pface_cells)
      if isone(length(lcells))
        pface_to_isboundary[k] = true
        pface_to_lcell[k] = Int8(lcells[1])
      end
      k += 1
    end
  end

  return pface_to_isboundary, pface_to_lcell
end

function get_pface_to_patch(ptopo::PatchTopology,Df::Integer)
  patch_faces = get_patch_faces(ptopo,Df)
  return Arrays.block_identity_array(Int32, patch_faces.ptrs)
end

function get_pface_to_lpface(ptopo::PatchTopology,Df::Integer)
  patch_faces = get_patch_faces(ptopo,Df)
  return Arrays.local_identity_array(Int32, patch_faces.ptrs)
end

function get_patch_to_tfaces(ptopo::PatchTopology,d::Integer,::IdentityVector)
  patch_to_faces = get_patch_faces(ptopo,d)
  Table(IdentityVector(num_faces(ptopo,d)),patch_to_faces.ptrs)
end

function get_patch_to_tfaces(ptopo::PatchTopology,d::Integer,tface_to_pface)
  pface_to_tface = find_inverse_index_map(tface_to_pface,num_faces(ptopo,d))
  patch_to_faces = get_patch_faces(ptopo,d)

  n_tfaces = length(tface_to_pface)
  n_patches = num_patches(ptopo)
  data = zeros(Int32, n_tfaces)
  ptrs = zeros(Int32, n_patches+1)

  k = 0
  for patch in 1:n_patches
    pfaces = patch_to_faces.ptrs[patch]:patch_to_faces.ptrs[patch+1]-1
    for pface in pfaces
      tface = pface_to_tface[pface]
      if tface > 0
        data[k+1] = tface
        ptrs[patch+1] += 1
        k += 1
      end
    end
  end
  Arrays.length_to_ptrs!(ptrs)

  return Table(data,ptrs)
end

function Base.view(ptopo::PatchTopology, cell_to_parent_cell::AbstractArray{<:Integer})
  ptopo_view, _ = patch_topology_view(ptopo, cell_to_parent_cell)
  return ptopo_view
end

function patch_topology_view(ptopo::PatchTopology, ::IdentityVector)
  D = num_cell_dims(ptopo)
  d_to_dpface_to_parent_dpface = [
    IdentityVector(num_faces(ptopo,d)) for d in 0:D
  ]
  return ptopo, d_to_dpface_to_parent_dpface
end

function patch_topology_view(ptopo::PatchTopology, cell_to_parent_cell::AbstractArray{<:Integer})
  Dc = num_cell_dims(ptopo)
  n_patches = num_patches(ptopo)

  topo = ptopo.topo
  cell_to_tag = fill(1,num_cells(topo))
  cell_to_tag[cell_to_parent_cell] .= 2
  labels = face_labeling_from_cell_tags(topo,cell_to_tag,["out","in"])

  d_to_patch_to_dfaces = Vector{Table{Int32,Vector{Int32},Vector{Int32}}}(undef,Dc+1)
  d_to_dpface_to_parent_dpface = Vector{Vector{Int32}}(undef,Dc+1)
  for d in 0:Dc
    patch_to_dfaces = get_patch_faces(ptopo,d)
    dface_mask = get_face_mask(labels,"in",d)
    is_in = Base.Fix1(getindex, dface_mask)

    ptrs = zeros(Int32, n_patches+1)
    for patch in 1:n_patches
      dfaces = view(patch_to_dfaces,patch)
      ptrs[patch+1] = count(is_in,dfaces)
    end
    Arrays.length_to_ptrs!(ptrs)

    data = zeros(Int32, ptrs[end]-1)
    for patch in 1:n_patches
      dfaces = view(patch_to_dfaces,patch)
      range = ptrs[patch]:ptrs[patch+1]-1
      data[range] .= Iterators.filter(is_in, dfaces)
    end

    d_to_patch_to_dfaces[d+1] = Table(data,ptrs)
    d_to_dpface_to_parent_dpface[d+1] = findall(lazy_map(Reindex(dface_mask),patch_to_dfaces.data))
  end

  ptopo_view = PatchTopology(topo, d_to_patch_to_dfaces, ptopo.metadata)
  return ptopo_view, d_to_dpface_to_parent_dpface
end

function restrict(ptopo::PatchTopology, ::IdentityVector; kwargs...)
  D = num_cell_dims(ptopo)
  d_to_dpface_to_parent_dpface = [
    IdentityVector(num_faces(ptopo,d)) for d in 0:D
  ]
  return ptopo, d_to_dpface_to_parent_dpface
end

function restrict(
  ptopo::PatchTopology, cell_to_parent_cell::AbstractArray{<:Integer}; 
  remove_empty_patches = false
)
  Dc = num_cell_dims(ptopo)
  n_patches = num_patches(ptopo)

  topo = ptopo.topo
  restricted_topo, d_to_dface_to_parent_dface = restrict(topo,cell_to_parent_cell)

  d_to_patch_to_dfaces = Vector{Table{Int32,Vector{Int32},Vector{Int32}}}(undef,Dc+1)
  for d in 0:Dc
    dface_to_parent_dface = d_to_dface_to_parent_dface[d+1]
    parent_dface_to_dface = find_inverse_index_map(
      dface_to_parent_dface, num_faces(topo,d)
    )
    patch_to_parent_dfaces = get_patch_faces(ptopo,d)
    
    ptrs = zeros(Int32, n_patches+1)
    for patch in 1:n_patches
      parent_dfaces = view(patch_to_parent_dfaces,patch)
      dfaces = view(parent_dface_to_dface,parent_dfaces)
      ptrs[patch+1] = count(!iszero,dfaces)
    end
    Arrays.length_to_ptrs!(ptrs)

    data = zeros(Int32, ptrs[end]-1)
    for patch in 1:n_patches
      parent_dfaces = view(patch_to_parent_dfaces,patch)
      dfaces = view(parent_dface_to_dface,parent_dfaces)
      range = ptrs[patch]:ptrs[patch+1]-1
      data[range] .= Iterators.filter(!iszero, dfaces)
    end

    d_to_patch_to_dfaces[d+1] = Table(data,ptrs)
    if remove_empty_patches
      Arrays.remove_empty_entries!(d_to_patch_to_dfaces[d+1])
    end
  end

  if isa(ptopo.metadata,StarPatchMetadata)
    root_dim = ptopo.metadata.root_dim
    parent_rface_to_rface = find_inverse_index_map(
      d_to_dface_to_parent_dface[root_dim+1],num_faces(topo,root_dim)
    )
    patch_roots = parent_rface_to_rface[ptopo.metadata.patch_roots]
    if remove_empty_patches
      patch_roots = filter!(iszero, patch_roots)
    end
    metadata = StarPatchMetadata(root_dim, patch_roots)
  else
    metadata = nothing
  end

  new_ptopo = PatchTopology(
    restricted_topo, d_to_patch_to_dfaces, metadata
  )
  return new_ptopo, d_to_dface_to_parent_dface
end

function extend_patches_by_single_layer(ptopo::PatchTopology{Dc}) where Dc
  topo = ptopo.topo
  nodes_to_cells = get_faces(topo,0,Dc)
  patch_to_nodes = get_patch_faces(ptopo,0)

  patch_cells = Arrays.merge_entries(
    nodes_to_cells, patch_to_nodes; 
    acc  = SortedSet{Int32}(), 
    post = dofs -> filter(x -> x > 0, dofs)
  )
  return PatchTopology(topo, patch_cells, nothing)
end

# PatchTriangulation

struct PatchGlue{Dc,A}
  tface_to_pface  :: A
  tface_to_patch  :: Vector{Int32}
  patch_to_tfaces :: Table{Int32,Vector{Int32},Vector{Int32}}
  function PatchGlue{Dc}(
    ptopo::PatchTopology,tface_to_pface
  ) where Dc
    pface_to_patch = get_pface_to_patch(ptopo,Dc)
    tface_to_patch = pface_to_patch[tface_to_pface]
    patch_to_tfaces = get_patch_to_tfaces(ptopo,Dc,tface_to_pface)

    A = typeof(tface_to_pface)
    new{Dc,A}(tface_to_pface,tface_to_patch,patch_to_tfaces)
  end
end

"""
    struct PatchTriangulation{Dc,Dp} <: Triangulation{Dc,Dp}
      trian :: Triangulation{Dc,Dp}
      ptopo :: PatchTopology
      glue  :: PatchGlue{Dc}
    end

Wrapper around a Triangulation, for patch-based assembly.

# Fields:
- `trian`: Underlying triangulation. In general, this can be a non-injective triangulation.
- `ptopo`: Patch topology, to which the triangulation faces can be mapped
- `glue`: Patch glue, mapping triangulation faces to the patches
"""
struct PatchTriangulation{Dc,Dp,A,B,C} <: Triangulation{Dc,Dp}
  trian :: A
  ptopo :: B
  glue  :: C
  function PatchTriangulation(
    trian::Triangulation{Dc,Dp},
    ptopo::PatchTopology,
    tface_to_pface::AbstractVector{<:Integer}
  ) where {Dc,Dp}
    glue = PatchGlue{Dc}(ptopo,tface_to_pface)
    A, B, C = typeof(trian), typeof(ptopo), typeof(glue)
    new{Dc,Dp,A,B,C}(trian,ptopo,glue)
  end
end

# Triangulation API

get_background_model(t::PatchTriangulation) = get_background_model(t.trian)
get_grid(t::PatchTriangulation) = get_grid(t.trian)
get_glue(t::PatchTriangulation,::Val{D}) where D = get_glue(t.trian,Val(D))
get_facet_normal(trian::PatchTriangulation) = get_facet_normal(trian.trian)

# Constructors

"""
    PatchTriangulation(model::DiscreteModel,ptopo::PatchTopology;tags=nothing)
    PatchTriangulation(::Type{ReferenceFE{D}},model::DiscreteModel,ptopo::PatchTopology;tags=nothing)
"""
function PatchTriangulation(model::DiscreteModel,ptopo::PatchTopology;kwargs...)
  D = min(num_cell_dims(model), num_cell_dims(ptopo))
  PatchTriangulation(ReferenceFE{D},model,ptopo;kwargs...)
end

function PatchTriangulation(
  ::Type{ReferenceFE{D}},model::DiscreteModel,ptopo::PatchTopology;tags=nothing
) where D
  @check D <= num_cell_dims(ptopo) <= num_cell_dims(model)
  patch_faces = get_patch_faces(ptopo,D)
  pface_to_face = patch_faces.data

  if !isnothing(tags)
    labels = get_face_labeling(model)
    face_to_mask = get_face_mask(labels,tags,D)
    pface_to_mask = lazy_map(Reindex(face_to_mask),pface_to_face)
    tface_to_pface = findall(pface_to_mask)
  else
    tface_to_pface = IdentityVector(num_faces(ptopo,D))
  end

  tface_to_face = lazy_map(Reindex(pface_to_face),tface_to_pface)
  trian = Triangulation(ReferenceFE{D},model,tface_to_face)
  return PatchTriangulation(trian,ptopo,tface_to_pface)
end

"""
    PatchBoundaryTriangulation(model::DiscreteModel,ptopo::PatchTopology{Dc};tags=nothing)
"""
function PatchBoundaryTriangulation(model::DiscreteModel{Dc},ptopo::PatchTopology{Dc};tags=nothing) where Dc
  patch_faces = get_patch_faces(ptopo,Dc-1)
  pface_to_face = patch_faces.data

  # Create p-face mask
  if !isnothing(tags)
    labels = get_face_labeling(model)
    face_to_mask = get_face_mask(labels,tags,Dc-1)
    pface_to_mask = lazy_map(Reindex(face_to_mask),pface_to_face)
  else
    pface_to_mask = Fill(true,num_faces(ptopo,Dc-1))
  end

  # Merge with boundary mask
  pface_to_isboundary, pface_to_lcell = get_patch_boundary_info(ptopo)
  tface_to_pface = findall(lazy_map(&,pface_to_isboundary,pface_to_mask))
  tface_to_lcell = pface_to_lcell[tface_to_pface]
  tface_to_face  = pface_to_face[tface_to_pface]

  # Create overlapping boundary triangulation
  topo = get_grid_topology(model)
  cell_grid = get_grid(model)
  face_grid = view(Grid(ReferenceFE{Dc-1},model),tface_to_face)
  face_glue = OverlappingFaceToCellGlue(topo,cell_grid,face_grid,tface_to_face,tface_to_lcell)
  face_trian = BodyFittedTriangulation(model,face_grid,tface_to_face)
  trian = BoundaryTriangulation(face_trian,face_glue)

  return PatchTriangulation(trian,ptopo,tface_to_pface)
end

"""
    PatchSkeletonTriangulation(model::DiscreteModel,ptopo::PatchTopology{Dc};tags=nothing)
"""
function PatchSkeletonTriangulation(model::DiscreteModel{Dc},ptopo::PatchTopology{Dc}; tags=nothing) where Dc
  patch_faces = get_patch_faces(ptopo,Dc-1)
  pface_to_face = patch_faces.data

  # Create p-face mask
  if !isnothing(tags)
    labels = get_face_labeling(model)
    face_to_mask = get_face_mask(labels,tags,Dc-1)
    pface_to_mask = lazy_map(Reindex(face_to_mask),pface_to_face)
  else
    pface_to_mask = Fill(true,num_faces(ptopo,Dc-1))
  end

  # Merge with boundary mask
  pface_to_isboundary, _ = get_patch_boundary_info(ptopo)
  tface_to_pface = findall(lazy_map((a,b) -> !a && b, pface_to_isboundary, pface_to_mask))
  tface_to_face  = pface_to_face[tface_to_pface]
  ntfaces = length(tface_to_face)

  # Create overlapping boundary triangulation
  topo = get_grid_topology(model)
  cell_grid = get_grid(model)
  face_grid = view(Grid(ReferenceFE{Dc-1},model),tface_to_face)
  face_trian = BodyFittedTriangulation(model,face_grid,tface_to_face)

  face_glue_plus = OverlappingFaceToCellGlue(topo,cell_grid,face_grid,tface_to_face,fill(Int8(1),ntfaces))
  trian_plus = BoundaryTriangulation(face_trian,face_glue_plus)

  face_glue_minus = OverlappingFaceToCellGlue(topo,cell_grid,face_grid,tface_to_face,fill(Int8(2),ntfaces))
  trian_minus = BoundaryTriangulation(face_trian,face_glue_minus)

  trian = SkeletonTriangulation(trian_plus,trian_minus)
  return PatchTriangulation(trian,ptopo,tface_to_pface)
end

for TT in (:PatchTriangulation,:PatchBoundaryTriangulation,:PatchSkeletonTriangulation)
  @eval begin
    function $TT(trian::Triangulation{Dc},ptopo::PatchTopology{Dc}; kwargs...) where Dc
      # Create the Boundary/Skeleton triangulation on a restricted PatchTopology
      model = get_background_model(trian)
      tcell_to_mcell = get_glue(trian,Val(Dc)).tface_to_mface
      ptopo_view, d_to_dface_to_parent_dface = patch_topology_view(ptopo,tcell_to_mcell)
      ptrian = $TT(model,ptopo_view; kwargs...)

      # Reindex the tface_to_pface map
      tface_to_pface = ptrian.glue.tface_to_pface
      tface_to_parent_pface = d_to_dface_to_parent_dface[Dc][tface_to_pface]
      return PatchTriangulation(
        ptrian.trian, ptopo, tface_to_parent_pface
      )
    end
  end
end

function get_patch_faces(trian::PatchTriangulation)
  Df = num_cell_dims(trian)
  tface_to_face = get_glue(trian,Val(Df)).tface_to_mface
  patch_to_tfaces = trian.glue.patch_to_tfaces
  return Table(lazy_map(Reindex(tface_to_face),patch_to_tfaces.data),patch_to_tfaces.ptrs)
end

# Domain changes between PatchTriangulations: 
#
# The idea here is that we cannot go from a PatchTriangulation to a Triangulation because
# the PatchTriangulation, in general, IS NOT injective. 
# I.e the tface_to_face map is not injective.
# However, two PatchTriangulations from the same PatchTopology can be mapped to each other
# since the map tface_to_pface IS injective.

function get_patch_glue(trian::PatchTriangulation)
  ptopo, glue = trian.ptopo, trian.glue

  tface_to_pface = glue.tface_to_pface
  pface_to_tface = PosNegPartition(tface_to_pface,num_faces(ptopo,num_cell_dims(trian)))
  tface_to_pface_map = Fill(GenericField(identity),num_cells(trian))
  FaceToFaceGlue(tface_to_pface,tface_to_pface_map,pface_to_tface)
end

function is_change_possible(strian::PatchTriangulation,ttrian::PatchTriangulation)
  if strian === ttrian || is_change_possible(strian.trian,ttrian.trian)
    return true
  end
  if num_cell_dims(strian) != num_cell_dims(ttrian)
    return false
  end
  sglue = get_patch_glue(strian)
  tglue = get_patch_glue(ttrian)
  return is_change_possible(sglue,tglue)
end
