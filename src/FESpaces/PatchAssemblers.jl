
using Gridap.Geometry: PatchTopology, get_patch_faces
using DataStructures: SortedSet
using BlockArrays

# PatchAssembler

struct PatchAssemblyStrategy{A,B} <: AssemblyStrategy
  ptopo :: PatchTopology
  patch_rows :: A
  patch_cols :: B
end

function map_patch_rows!(gids,a::PatchAssemblyStrategy,rows,patch)
  prows = a.patch_rows[patch]
  u = -one(eltype(gids))
  for i in eachindex(rows)
    ri = rows[i]
    if ri > 0
      gids[i] = searchsortedfirst(prows,ri)
    else
      gids[i] = u
    end
  end
  nothing
end

function map_patch_cols!(gids,a::PatchAssemblyStrategy,cols,patch)
  pcols = a.patch_cols[patch]
  u = -one(eltype(gids))
  for i in eachindex(cols)
    ci = cols[i]
    if ci > 0
      gids[i] = searchsortedfirst(pcols,ci)
    else
      gids[i] = u
    end
  end
  nothing
end

Arrays.return_cache(a::AssemblyStrategyMap,ids::AbstractArray,patch) = return_cache(a,ids)

function Arrays.evaluate!(cache,a::AssemblyStrategyMap{:cols},ids::AbstractArray,patch)
  setsize!(cache,size(ids))
  gids = cache.array
  map_patch_cols!(gids,a.strategy,ids,patch)
  gids
end

function Arrays.evaluate!(cache,a::AssemblyStrategyMap{:rows},ids::AbstractArray,patch)
  setsize!(cache,size(ids))
  gids = cache.array
  map_patch_rows!(gids,a.strategy,ids,patch)
  gids
end

Arrays.return_cache(k::AssemblyStrategyMap,ids::ArrayBlock,patch) = return_cache(k,ids)

function Arrays.evaluate!(cache,k::AssemblyStrategyMap,ids::ArrayBlock,patch)
  a,b = cache
  for i in eachindex(ids.array)
    if ids.touched[i]
      a.array[i] = evaluate!(b[i],k,ids.array[i],patch)
    end
  end
  a
end

struct PatchAssembler <: Assembler
  ptopo    :: PatchTopology
  strategy :: PatchAssemblyStrategy
  rows
  cols
end

get_rows(assembler::PatchAssembler) = assembler.rows
get_cols(assembler::PatchAssembler) = assembler.cols

function PatchAssembler(ptopo::PatchTopology,trial::FESpace,test::FESpace)
  patch_rows = get_patch_dofs(test,ptopo)
  patch_cols = get_patch_dofs(trial,ptopo)
  strategy = PatchAssemblyStrategy(ptopo,patch_rows,patch_cols)
  rows = blockedrange(map(length,patch_rows))
  cols = blockedrange(map(length,patch_cols))
  return PatchAssembler(ptopo,strategy,rows,cols)
end

function get_patch_dofs(space::FESpace,ptopo::PatchTopology)
  trian = get_triangulation(space)

  Df = num_cell_dims(trian)
  face_to_tface = get_glue(trian,Val(Df)).mface_to_tface
  @notimplementedif isnothing(face_to_tface)
  tface_to_dofs = get_cell_dof_ids(space)
  patch_to_faces = get_patch_faces(ptopo,Df)

  patch_to_dofs = map(patch_to_faces) do pfaces
    tfaces = filter(!iszero,face_to_tface[pfaces])
    dofs = SortedSet{Int}()
    for tface in tfaces
      push!(dofs,view(tface_to_dofs,tface)...)
    end
    collect(dofs)
  end |> Table

  return patch_to_dofs
end

"""
function get_cell_dof_ids(f::FESpace,ttrian::PatchTriangulation)
  tface_to_dofs = get_cell_dof_ids(f,ttrian.trian)

  trian = get_triangulation(f)
  pface_to_patch = get_pface_to_patch(ttrian.ptopo,num_cell_dims(trian))
  tface_to_patch = view(pface_to_patch,ttrian.tface_to_pface)

  return pair_arrays(tface_to_dofs,tface_to_patch)
end
"""

"""
function get_patch_dofs(space::MultiFieldFESpace,ptopo::PatchTopology)
  offsets = get_field_offsets(space)
  sf_patch_dofs = map(space) do sf
    get_patch_dofs(sf,ptopo)
  end
  mf_patch_dofs = append_tables_locally(offsets,sf_patch_dofs)
  return mf_patch_dofs
end
"""

function assemble_matrix(assem::PatchAssembler,cellmat)
  patch_sizes = map(tuple,blocklengths(assem.rows),blocklengths(assem.cols))
  PatchAssemblyMap(assem,patch_sizes,cellmat)
end

function assemble_vector(assem::PatchAssembler,cellvec)
  patch_sizes = map(tuple,blocklengths(assem.rows))
  PatchAssemblyMap(assem,patch_sizes,cellvec)
end

function assemble_matrix_and_vector(assem::PatchAssembler,celldata)
  patch_sizes = map((r,c) -> ((r,c),r),blocklengths(assem.rows),blocklengths(assem.cols))
  PatchAssemblyMap(assem,patch_sizes,celldata)
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

function collect_cell_patch(ptopo::PatchTopology,a::DomainContribution)
  p, q = [], []
  for strian in get_domains(a)
    @assert isa(strian,Geometry.PatchTriangulation)
    @assert strian.ptopo === ptopo
    Df = num_cell_dims(strian)

    patch_to_faces = get_patch_faces(ptopo,Df)
    pface_to_tface = find_inverse_index_map(strian.tface_to_pface,num_faces(ptopo,Df))
    patch_to_tfaces = Table(pface_to_tface,patch_to_faces.ptrs)
    push!(p,patch_to_tfaces)

    pface_to_patch = Geometry.get_pface_to_patch(ptopo,Df)
    tface_to_patch = lazy_map(Reindex(pface_to_patch),strian.tface_to_pface)
    push!(q,tface_to_patch)
  end
  return p, q
end

function attach_patch_map(k,r,q)
  c = []
  for (ri,qi) in zip(r,q)
    push!(c,lazy_map(k,ri,qi))
  end
  return c
end

function collect_patch_cell_matrix(assem::PatchAssembler,trial::FESpace,test::FESpace,a::DomainContribution)
  w, _r, _c = collect_cell_matrix(trial,test,a)
  p, q = collect_cell_patch(assem.ptopo,a)
  r = attach_patch_map(AssemblyStrategyMap{:rows}(assem.strategy),_r,q)
  c = attach_patch_map(AssemblyStrategyMap{:cols}(assem.strategy),_c,q)
  return (w, r, c, p)
end

function collect_patch_cell_vector(assem::PatchAssembler,test::FESpace,a::DomainContribution)
  w, r = collect_cell_vector(test,a)
  p, q = collect_cell_patch(assem.ptopo,a)
  r = attach_patch_map(AssemblyStrategyMap{:rows}(assem.strategy),_r,q)
  return (w, r, p)
end

function _collect_patch_cell_matvec(assem::PatchAssembler,trial::FESpace,test::FESpace,a::DomainContribution)
  w, r, c = _collect_cell_matvec(trial,test,a)
  p, q = collect_cell_patch(assem.ptopo,a)
  r = attach_patch_map(AssemblyStrategyMap{:rows}(assem.strategy),_r,q)
  c = attach_patch_map(AssemblyStrategyMap{:cols}(assem.strategy),_c,q)
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

# PatchAssemblyMap

struct PatchAssemblyMap{T} <: Map
  assem :: PatchAssembler
  patch_sizes :: Vector{T}
  cell_data
end

function Arrays.return_cache(k::PatchAssemblyMap{NTuple{2,Int}},patch)
  res = zeros(k.patch_sizes[patch])

  caches = ()
  for (cell_vals, cell_rows, cell_cols, patch_cells) in zip(k.cell_data...)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    vals_cache = array_cache(cell_vals)

    rows1 = getindex!(rows_cache,cell_rows,1)
    cols1 = getindex!(cols_cache,cell_cols,1)
    vals1 = getindex!(vals_cache,cell_vals,1)

    add_cache = return_cache(AddEntriesMap(+),res,vals1,rows1,cols1)

    add_caches = (add_cache, vals_cache, rows_cache, cols_cache)
    cells_cache = array_cache(patch_cells)
    caches = (caches...,(add_caches, cells_cache))
  end

  return CachedArray(res), caches
end

function Arrays.evaluate!(cache,k::PatchAssemblyMap{NTuple{2,Int}},patch)
  c_res, caches = cache

  setsize!(c_res,k.patch_sizes[patch])
  res = c_res.array

  fill!(res,zero(eltype(res)))
  
  for (c,cell_vals, cell_rows, cell_cols, patch_cells) in zip(caches, k.cell_data...)
    ac, cc = c
    cells = getindex!(cc, patch_cells, patch)
    _numeric_loop_matrix!(res,ac,cell_vals,cell_rows,cell_cols,cells)
  end
  
  return res
end
