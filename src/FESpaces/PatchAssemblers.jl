
using Gridap.Geometry: PatchTopology, get_patch_faces
using DataStructures: SortedSet
using BlockArrays

# PatchAssembler
struct PatchAssemblyStrategy{A,B} <: AssemblyStrategy
  ptopo :: PatchTopology
  patch_rows :: A
  patch_cols :: B
end

function map_rows!(gids,a::PatchAssemblyStrategy,rows,patch)
  prows = view(a.patch_rows,patch)
  n = length(prows)+1
  u = -one(eltype(gids))
  for i in eachindex(rows)
    ri = rows[i]
    if ri > 0
      gid = searchsortedfirst(prows,ri)
      found = !isequal(gid,n) && isequal(prows[gid],ri)
      gids[i] = ifelse(found, gid, u)
    else
      gids[i] = u
    end
  end
  nothing
end

function map_cols!(gids,a::PatchAssemblyStrategy,cols,patch)
  pcols = view(a.patch_cols,patch)
  n = length(pcols)+1
  u = -one(eltype(gids))
  for i in eachindex(cols)
    ci = cols[i]
    if ci > 0
      gid = searchsortedfirst(pcols,ci)
      found = !isequal(gid,n) && isequal(pcols[gid],ci)
      gids[i] = ifelse(found, gid, u)
    else
      gids[i] = u
    end
  end
  nothing
end

struct PatchAssembler <: Assembler
  ptopo :: PatchTopology
  strategy
  rows
  cols
end

get_rows(assembler::PatchAssembler) = assembler.rows
get_cols(assembler::PatchAssembler) = assembler.cols

function PatchAssembler(ptopo::PatchTopology,trial::FESpace,test::FESpace;kwargs...)
  patch_rows = get_patch_assembly_ids(test,ptopo;kwargs...)
  patch_cols = get_patch_assembly_ids(trial,ptopo;kwargs...)
  strategy = PatchAssemblyStrategy(ptopo,patch_rows,patch_cols)
  rows = map(length,patch_rows)
  cols = map(length,patch_cols)
  return PatchAssembler(ptopo,strategy,rows,cols)
end

function get_patch_assembly_ids(space::FESpace,ptopo::PatchTopology;assembly=:all)
  if assembly == :all
    _patch_assembly_ids_all(space,ptopo)
  elseif assembly == :star
    # Takes all dofs owned by a face that is connected to the root of the patch.
    # This is different than :interior for patches whose root is on the domain boundary.
    _patch_assembly_ids_star(space,ptopo)
  elseif assembly == :interior
    # Takes all dofs owned by a face that is not boundary of the patch
    _patch_assembly_ids_interior(space,ptopo)
  elseif assembly == :boundary
    # Takes all dofs owned by a face that is boundary of the patch
    _patch_assembly_ids_boundary(space,ptopo)
  else
    @notimplemented """
    Assembly type $assembly not implemented. Options are (:all, :star, :interior, :boundary).
    """
  end :: Table{Int32,Vector{Int32},Vector{Int32}}
end

function _patch_assembly_ids_all(space::FESpace,ptopo::PatchTopology)
  trian = get_triangulation(space)
  Df = num_cell_dims(trian)
  face_to_tface = get_glue(trian,Val(Df)).mface_to_tface
  @notimplementedif isnothing(face_to_tface)

  face_dof_ids = extend(get_cell_dof_ids(space),face_to_tface)
  patch_to_faces = get_patch_faces(ptopo,Df)

  patch_rows = Arrays.merge_entries(
    face_dof_ids, patch_to_faces ; 
    acc  = SortedSet{Int32}(), 
    post = dofs->filter(x -> x > 0, dofs)
  )
  return patch_rows
end

function _patch_assembly_ids_star(space::FESpace, ptopo::PatchTopology)
  @check ptopo.metadata isa Geometry.StarPatchMetadata """
    PatchTopology does not have StarPatchMetadata metadata.
  """
  D = num_cell_dims(ptopo)
  topo, metadata = ptopo.topo, ptopo.metadata
  Dr, patch_roots = metadata.root_dim, metadata.patch_roots

  # A dpface is masked (removed) iff
  #  - It is a boundary of the patch, AND
  #  - It is not connected to the root of the patch
  d_to_dpface_to_mask = map(0:D) do d

    dpface_to_dface = Geometry.get_patch_faces(ptopo,d).data
    dpface_to_patch = Geometry.get_pface_to_patch(ptopo,d)
    root_to_dfaces = Geometry.get_faces(topo,Dr,d)
    
    dpface_to_mask = Geometry.compute_isboundary_face(ptopo,d)
    for (dpface, dface) in enumerate(dpface_to_dface)
      !dpface_to_mask[dpface] && continue # Already interior
      # Check if it is connected to the root
      root = patch_roots[dpface_to_patch[dpface]]
      dpface_to_mask[dpface] = (dface ∉ view(root_to_dfaces,root))
    end

    return dpface_to_mask
  end

  return _patch_assembly_ids_masked(space,ptopo,d_to_dpface_to_mask)
end

function _patch_assembly_ids_boundary(space::FESpace, ptopo::PatchTopology;reverse=false)
  # A dpface is masked (removed) iff it is NOT boundary of the patch
  D = num_cell_dims(ptopo)
  d_to_dpface_to_mask = map(0:D) do d
    Geometry.compute_isboundary_face(ptopo,d)
  end
  return _patch_assembly_ids_masked(space,ptopo,d_to_dpface_to_mask;reverse)
end

function _patch_assembly_ids_interior(space::FESpace, ptopo::PatchTopology)
  return _patch_assembly_ids_boundary(space,ptopo;reverse=true)
end

function _patch_assembly_ids_masked(
  space::FESpace, ptopo::PatchTopology, d_to_dpface_to_mask;
  reverse = true
)
  strian = get_triangulation(space)
  model = get_background_model(strian)

  D = num_cell_dims(strian)
  ttrian = PatchTriangulation(ReferenceFE{D},model,ptopo)

  tface_to_mface = get_glue(ttrian, Val(D)).tface_to_mface
  mface_to_sface = get_glue(strian, Val(D)).mface_to_tface
  tface_to_sface = mface_to_sface[tface_to_mface]

  cell_conformity = get_cell_conformity(space)
  d_to_tcell_to_tdface = [ Geometry.generate_patch_faces(ptopo,D,d) for d in 0:D ]
  tcell_to_ldof_mask = generate_cell_dof_mask(
    cell_conformity,tface_to_sface,d_to_tcell_to_tdface,d_to_dpface_to_mask;reverse
  )

  tcell_dof_ids = get_cell_dof_ids(space,ttrian)
  tcell_dof_ids_masked = lazy_map(getindex, tcell_dof_ids, tcell_to_ldof_mask)
  patch_to_tfaces = Geometry.get_patch_to_tfaces(ptopo,D,ttrian.glue.tface_to_pface)

  patch_rows = Arrays.merge_entries(
    tcell_dof_ids_masked, patch_to_tfaces; 
    acc  = SortedSet{Int32}(), 
    post = dofs -> filter(x -> x > 0, dofs)
  )

  return patch_rows
end

function assemble_matrix(assem::PatchAssembler,cellmat)
  patch_sizes = map(tuple,assem.rows,assem.cols)
  k = PatchAssemblyMap(assem,patch_sizes,cellmat)
  return lazy_map(k,collect(1:Geometry.num_patches(assem.ptopo)))
end

function assemble_vector(assem::PatchAssembler,cellvec)
  patch_sizes = map(tuple,assem.rows)
  k = PatchAssemblyMap(assem,patch_sizes,cellvec)
  return lazy_map(k,collect(1:Geometry.num_patches(assem.ptopo)))
end

function assemble_matrix_and_vector(assem::PatchAssembler,celldata)
  patch_sizes = map((r,c) -> ((r,c),(r,)),assem.rows,assem.cols)
  k = PatchAssemblyMap(assem,patch_sizes,celldata)
  return lazy_map(k,collect(1:Geometry.num_patches(assem.ptopo)))
end

function assemble_matrix(f::Function,a::PatchAssembler,U::FESpace,V::FESpace)
  v = get_fe_basis(V)
  u = get_trial_fe_basis(U)
  assemble_matrix(a,collect_patch_cell_matrix(a,U,V,f(u,v)))
end

function assemble_vector(f::Function,a::PatchAssembler,V::FESpace)
  v = get_fe_basis(V)
  assemble_vector(a,collect_patch_cell_vector(a,V,f(v)))
end

function assemble_matrix_and_vector(f::Function,b::Function,a::PatchAssembler,U::FESpace,V::FESpace)
  v = get_fe_basis(V)
  u = get_trial_fe_basis(U)
  assemble_matrix_and_vector(a,collect_patch_cell_matrix_and_vector(a,U,V,f(u,v),b(v)))
end

function assemble_matrix_and_vector(f::Function,b::Function,a::PatchAssembler,U::FESpace,V::FESpace,uhd)
  v = get_fe_basis(V)
  u = get_trial_fe_basis(U)
  assemble_matrix_and_vector(a,collect_patch_cell_matrix_and_vector(a,U,V,f(u,v),b(v),uhd))
end

function collect_cell_patch(ptopo::PatchTopology,a::DomainContribution)
  p, q = [], []
  for strian in get_domains(a)
    @assert isa(strian,Geometry.PatchTriangulation)
    @assert strian.ptopo === ptopo
    glue = strian.glue
    push!(p,glue.patch_to_tfaces)
    push!(q,glue.tface_to_patch)
  end
  return p, q
end

function attach_patch_map(f,strat,r,q)
  c = []
  for (ri,qi) in zip(r,q)
    push!(c,f(strat,ri,qi))
  end
  return c
end

function collect_patch_cell_matrix(assem::PatchAssembler,trial::FESpace,test::FESpace,a::DomainContribution)
  w, _r, _c = collect_cell_matrix(trial,test,a)
  p, q = collect_cell_patch(assem.ptopo,a)
  r = attach_patch_map(map_cell_rows,assem.strategy,_r,q)
  c = attach_patch_map(map_cell_cols,assem.strategy,_c,q)
  return (w, r, c, p)
end

function collect_patch_cell_vector(assem::PatchAssembler,test::FESpace,a::DomainContribution)
  w, _r = collect_cell_vector(test,a)
  p, q = collect_cell_patch(assem.ptopo,a)
  r = attach_patch_map(map_cell_rows,assem.strategy,_r,q)
  return (w, r, p)
end

function _collect_patch_cell_matvec(assem::PatchAssembler,trial::FESpace,test::FESpace,a::DomainContribution)
  w, _r, _c = _collect_cell_matvec(trial,test,a)
  p, q = collect_cell_patch(assem.ptopo,a)
  r = attach_patch_map(map_cell_rows,assem.strategy,_r,q)
  c = attach_patch_map(map_cell_cols,assem.strategy,_c,q)
  return (w, r, c, p)
end

function collect_patch_cell_matrix_and_vector(
  assem::PatchAssembler,trial::FESpace,test::FESpace,
  biform::DomainContribution,liform::DomainContribution
)
  matvec, mat, vec = _pair_contribution_when_possible(biform,liform)
  matvecdata = _collect_patch_cell_matvec(assem,trial,test,matvec)
  matdata = collect_patch_cell_matrix(assem,trial,test,mat)
  vecdata = collect_patch_cell_vector(assem,test,vec)
  (matvecdata, matdata, vecdata)
end

function collect_patch_cell_matrix_and_vector(
  assem::PatchAssembler,trial::FESpace,test::FESpace,
  biform::DomainContribution,liform::DomainContribution,uhd::FEFunction
)
  matvec, mat, vec = _pair_contribution_when_possible(biform,liform,uhd)
  matvecdata = _collect_patch_cell_matvec(assem,trial,test,matvec)
  matdata = collect_patch_cell_matrix(assem,trial,test,mat)
  vecdata = collect_patch_cell_vector(assem,test,vec)
  (matvecdata, matdata, vecdata)
end

function collect_and_merge_cell_matrix(assem::PatchAssembler,contributions...)
  data = ()
  for c in contributions
    data = (data..., collect_patch_cell_matrix(assem,c...))
  end
  merge_assembly_data(data...)
end

function collect_and_merge_cell_vector(assem::PatchAssembler,contributions...)
  data = ()
  for c in contributions
    data = (data..., collect_patch_cell_vector(assem,c...))
  end
  merge_assembly_data(data...)
end

function collect_and_merge_cell_matrix_and_vector(assem::PatchAssembler,contributions...)
  data = ()
  for c in contributions
    data = (data..., collect_patch_cell_matrix_and_vector(assem,c...))
  end
  merge_assembly_matvec_data(data...)
end

# PatchAssemblyMap

struct PatchAssemblyMap{T} <: Map
  assem :: PatchAssembler
  patch_sizes :: Vector{T}
  cell_data
end

_unview(a) = a
_unview(a::ArrayBlockView) = a.array

function _alloc_cache(::PatchAssemblyStrategy,s::NTuple{N,Int}) where N
  zeros(s...)
end

function _alloc_cache(::PatchAssemblyStrategy,s::Tuple{Vararg{BlockedOneTo}})
  N  = length(s)
  bs = map(blocklength,s)
  ss = map(blocklengths,s)
  array = [ zeros(ntuple(i -> ss[i][I[i]], Val(N))) for I in CartesianIndices(bs) ]
  ArrayBlock(array,fill(true,bs))
end

function _alloc_cache(b::ArrayBlockView,s::Tuple{Vararg{BlockedOneTo}})
  N = length(s)
  bs = map(blocklength,s)
  ss = map(blocklengths,s)
  array = [ zeros(ntuple(i -> ss[i][I[i]], Val(N))) for I in CartesianIndices(bs) ]
  bmap = ifelse(N == 2, b.block_map, map(idx -> CartesianIndex(idx[1]), diag(b.block_map)))
  ArrayBlockView(ArrayBlock(array,fill(true,bs)),bmap)
end

function _resize_cache!(a,::PatchAssemblyStrategy,s::NTuple{N,Int}) where N
  setsize!(a,s)
end

function _resize_cache!(a,::PatchAssemblyStrategy,s::Tuple{Vararg{BlockedOneTo}})
  N = length(s)
  bs = map(blocklength,s)
  ss = map(blocklengths,s)
  for I in CartesianIndices(bs)
    setsize!(a.array[I],ntuple(i -> ss[i][I[i]], Val(N)))
  end
end

function _resize_cache!(a,::ArrayBlockView,s::Tuple{Vararg{BlockedOneTo}})
  N = length(s)
  bs = map(blocklength,s)
  ss = map(blocklengths,s)
  for I in CartesianIndices(bs)
    setsize!(a.array.array[I],ntuple(i -> ss[i][I[i]], Val(N)))
  end
end

# Mat & Vec assembly

function Arrays.return_cache(k::PatchAssemblyMap,patch)
  res = _alloc_cache(k.assem.strategy,k.patch_sizes[patch])
  caches = patch_assembly_cache(res,k.cell_data)

  c_res = CachedArray(res)
  uwr_cache = return_cache(Arrays.unwrap_cached_array,c_res)
  return c_res, uwr_cache, caches
end

function Arrays.evaluate!(cache,k::PatchAssemblyMap,patch)
  c_res, uwr_cache, caches = cache
  _resize_cache!(c_res,k.assem.strategy,k.patch_sizes[patch])
  res = evaluate!(uwr_cache,Arrays.unwrap_cached_array,c_res)
  Arrays.fill_entries!(res,zero(Arrays.entry_type(res)))
  patch_assembly!(caches,res,k.cell_data,patch)
  return _unview(res)
end

# Mat-Vec assembly

function Arrays.return_cache(k::PatchAssemblyMap{<:Tuple{<:Tuple,<:Tuple}},patch)
  mat = _alloc_cache(k.assem.strategy,k.patch_sizes[patch][1])
  vec = _alloc_cache(k.assem.strategy,k.patch_sizes[patch][2])

  matvecdata, matdata, vecdata = k.cell_data
  mat_caches = patch_assembly_cache(mat,matdata)
  vec_caches = patch_assembly_cache(vec,vecdata)
  matvec_caches = patch_assembly_cache(mat,vec,matvecdata)

  c_mat, c_vec = CachedArray(mat), CachedArray(vec)
  uwm_cache = return_cache(Arrays.unwrap_cached_array,c_mat)
  uwv_cache = return_cache(Arrays.unwrap_cached_array,c_vec)

  return c_mat, c_vec, uwm_cache, uwv_cache, matvec_caches, mat_caches, vec_caches
end

function Arrays.evaluate!(cache,k::PatchAssemblyMap{<:Tuple{<:Tuple,<:Tuple}},patch)
  c_mat, c_vec, uwm_cache, uwv_cache, matvec_caches, mat_caches, vec_caches = cache
  matvecdata, matdata, vecdata = k.cell_data

  _resize_cache!(c_mat,k.assem.strategy,k.patch_sizes[patch][1])
  _resize_cache!(c_vec,k.assem.strategy,k.patch_sizes[patch][2])
  mat = evaluate!(uwm_cache,Arrays.unwrap_cached_array,c_mat)
  vec = evaluate!(uwv_cache,Arrays.unwrap_cached_array,c_vec)

  Arrays.fill_entries!(mat,zero(Arrays.entry_type(mat)))
  Arrays.fill_entries!(vec,zero(Arrays.entry_type(vec)))

  patch_assembly!(matvec_caches,mat,vec,matvecdata,patch)
  patch_assembly!(mat_caches,mat,matdata,patch)
  patch_assembly!(vec_caches,vec,vecdata,patch)
  
  return _unview(mat), _unview(vec)
end

const MatOrMatBlock = Union{AbstractMatrix,MatrixBlock,MatrixBlockView}
const VecOrVecBlock = Union{AbstractVector,VectorBlock,VectorBlockView}

function patch_assembly_cache(mat::MatOrMatBlock,cell_matdata)
  caches = ()
  for (cell_vals, cell_rows, cell_cols, patch_cells) in zip(cell_matdata...)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    vals_cache = array_cache(cell_vals)

    rows1 = getindex!(rows_cache,cell_rows,1)
    cols1 = getindex!(cols_cache,cell_cols,1)
    vals1 = getindex!(vals_cache,cell_vals,1)

    add_cache = return_cache(AddEntriesMap(+),mat,vals1,rows1,cols1)

    add_caches = (add_cache, vals_cache, rows_cache, cols_cache)
    cells_cache = array_cache(patch_cells)
    caches = (caches...,(add_caches, cells_cache))
  end
  return caches
end

function patch_assembly!(caches,mat::MatOrMatBlock,cell_matdata,patch)
  for (c, cell_vals, cell_rows, cell_cols, patch_cells) in zip(caches, cell_matdata...)
    ac, cc = c
    cells = getindex!(cc, patch_cells, patch)
    _numeric_loop_matrix!(mat,ac,cell_vals,cell_rows,cell_cols,cells)
  end
end

function patch_assembly_cache(vec::VecOrVecBlock,cell_vecdata)
  caches = ()
  for (cell_vals, cell_rows, patch_cells) in zip(cell_vecdata...)
    rows_cache = array_cache(cell_rows)
    vals_cache = array_cache(cell_vals)

    rows1 = getindex!(rows_cache,cell_rows,1)
    vals1 = getindex!(vals_cache,cell_vals,1)

    add_cache = return_cache(AddEntriesMap(+),vec,vals1,rows1)

    add_caches = (add_cache, vals_cache, rows_cache)
    cells_cache = array_cache(patch_cells)
    caches = (caches...,(add_caches, cells_cache))
  end
  return caches
end

function patch_assembly!(caches,vec::VecOrVecBlock,cell_vecdata,patch)
  for (c, cell_vals, cell_rows, patch_cells) in zip(caches, cell_vecdata...)
    ac, cc = c
    cells = getindex!(cc, patch_cells, patch)
    _numeric_loop_vector!(vec,ac,cell_vals,cell_rows,cells)
  end
end

function patch_assembly_cache(mat::MatOrMatBlock,vec::VecOrVecBlock,cell_matvecdata)
  caches = ()
  for (cell_vals, cell_rows, cell_cols, patch_cells) in zip(cell_matvecdata...)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    vals_cache = array_cache(cell_vals)

    rows1 = getindex!(rows_cache,cell_rows,1)
    cols1 = getindex!(cols_cache,cell_cols,1)
    mat1, vec1 = getindex!(vals_cache,cell_vals,1)

    add_mat_cache = return_cache(AddEntriesMap(+),mat,mat1,rows1,cols1)
    add_vec_cache = return_cache(AddEntriesMap(+),vec,vec1,rows1)

    add_caches = (add_mat_cache, add_vec_cache, vals_cache, rows_cache, cols_cache)
    cells_cache = array_cache(patch_cells)
    caches = (caches...,(add_caches, cells_cache))
  end
  return caches
end

function patch_assembly!(caches,mat::MatOrMatBlock,vec::VecOrVecBlock,cell_matvecdata,patch)
  for (c, cell_vals, cell_rows, cell_cols, patch_cells) in zip(caches, cell_matvecdata...)
    ac, cc = c
    cells = getindex!(cc, patch_cells, patch)
    _numeric_loop_matvec!(mat,vec,ac,cell_vals,cell_rows,cell_cols,cells)
  end
end
