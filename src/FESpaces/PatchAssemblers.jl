
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

function map_cols!(gids,a::PatchAssemblyStrategy,cols,patch)
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

struct PatchAssembler <: Assembler
  ptopo :: PatchTopology
  strategy
  rows
  cols
end

get_rows(assembler::PatchAssembler) = assembler.rows
get_cols(assembler::PatchAssembler) = assembler.cols

function PatchAssembler(ptopo::PatchTopology,trial::FESpace,test::FESpace)
  patch_rows = get_patch_assembly_ids(test,ptopo)
  patch_cols = get_patch_assembly_ids(trial,ptopo)
  strategy = PatchAssemblyStrategy(ptopo,patch_rows,patch_cols)
  rows = map(length,patch_rows)
  cols = map(length,patch_cols)
  return PatchAssembler(ptopo,strategy,rows,cols)
end

function get_patch_assembly_ids(space::FESpace,ptopo::PatchTopology)
  trian = get_triangulation(space)
  Df = num_cell_dims(trian)
  face_to_tface = get_glue(trian,Val(Df)).mface_to_tface
  @notimplementedif isnothing(face_to_tface)

  face_dof_ids = extend(get_cell_dof_ids(space),face_to_tface)
  patch_to_faces = get_patch_faces(ptopo,Df)

  patch_rows = Arrays.merge_entries(
    face_dof_ids, patch_to_faces ; 
    acc  = SortedSet{Int32}(), 
    post = dofs->filter(x->x>0,dofs)
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
  uwr_cache = return_cache(Fields.unwrap_cached_array,c_res)
  return c_res, uwr_cache, caches
end

function Arrays.evaluate!(cache,k::PatchAssemblyMap,patch)
  c_res, uwr_cache, caches = cache
  _resize_cache!(c_res,k.assem.strategy,k.patch_sizes[patch])
  res = evaluate!(uwr_cache,Fields.unwrap_cached_array,c_res)
  Fields._zero_entries!(res)
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
  uwm_cache = return_cache(Fields.unwrap_cached_array,c_mat)
  uwv_cache = return_cache(Fields.unwrap_cached_array,c_vec)

  return c_mat, c_vec, uwm_cache, uwv_cache, matvec_caches, mat_caches, vec_caches
end

function Arrays.evaluate!(cache,k::PatchAssemblyMap{<:Tuple{<:Tuple,<:Tuple}},patch)
  c_mat, c_vec, uwm_cache, uwv_cache, matvec_caches, mat_caches, vec_caches = cache
  matvecdata, matdata, vecdata = k.cell_data

  _resize_cache!(c_mat,k.assem.strategy,k.patch_sizes[patch][1])
  _resize_cache!(c_vec,k.assem.strategy,k.patch_sizes[patch][2])
  mat = evaluate!(uwm_cache,Fields.unwrap_cached_array,c_mat)
  vec = evaluate!(uwv_cache,Fields.unwrap_cached_array,c_vec)

  Fields._zero_entries!(mat)
  Fields._zero_entries!(vec)

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

###########################################################################################

struct LocalOperator
  local_map :: Map
  trian_out :: Triangulation
  space_out :: FESpace
  weakform  :: Function
  collect_coefficients :: Bool
end

function LocalOperator(
  local_map::Map,
  ptopo::PatchTopology,
  space_to::FESpace,
  space_from::FESpace,
  lhs::Function,
  rhs::Function;
  space_out::FESpace = space_to,
  trian_out::Triangulation = get_triangulation(space_out),
  collect_coefficients::Bool = true
)
  function weakform(u_from)
    u_to = get_trial_fe_basis(space_to)
    v_to = get_fe_basis(space_to)

    lhs_assem = PatchAssembler(ptopo,space_to,space_to)
    lhs_mats = assemble_matrix(
      lhs_assem,collect_patch_cell_matrix(lhs_assem,space_to,space_to,lhs(u_to,v_to))
    )

    rhs_assem = PatchAssembler(ptopo,space_from,space_to)
    rhs_mats = assemble_matrix(
      rhs_assem,collect_patch_cell_matrix(rhs_assem,space_from,space_to,rhs(u_from,v_to))
    )

    pair_arrays(lhs_mats,rhs_mats)
  end

  return LocalOperator(local_map,trian_out,space_out,weakform,collect_coefficients)
end

function LocalOperator(
  local_map::Map,
  space_to::FESpace,
  lhs::Function,
  rhs::Function;
  space_out::FESpace = space_to,
  trian_out::Triangulation = get_triangulation(space_out),
  collect_coefficients::Bool = true
)
  function weakform(u_from)
    u_to = get_trial_fe_basis(space_to)
    v_to = get_fe_basis(space_to)

    lhs_c = lhs(u_to,v_to)
    @check all(t -> t === trian_out, get_domains(lhs_c))
    lhs_mats = get_contribution(lhs_c, trian_out)

    rhs_c = rhs(u_from,v_to)
    @check all(t -> t === trian_out, get_domains(rhs_c))
    rhs_mats = get_contribution(rhs_c, trian_out)

    pair_arrays(lhs_mats,rhs_mats)
  end

  return LocalOperator(local_map,trian_out,space_out,weakform,collect_coefficients)
end

(P::LocalOperator)(u) = evaluate(P,u)

function Arrays.evaluate!(
  cache,k::LocalOperator,space::FESpace
)
  u = evaluate!(cache,k,get_trial_fe_basis(space))
  v = similar_cell_field(u,lazy_map(transpose,CellData.get_data(u)),get_triangulation(u),DomainStyle(u))
  return u, v
end

function Arrays.evaluate!(
  cache,k::LocalOperator,v::SingleFieldFEBasis{<:TestBasis}
)
  u = FESpaces.similar_fe_basis(
    v,lazy_map(transpose,get_data(v)),get_triangulation(v),TrialBasis(),DomainStyle(v)
  )
  data = _compute_local_solves(k,u)
  return GenericCellField(data,k.trian_out,ReferenceDomain())
end

function Arrays.evaluate!(
  cache,k::LocalOperator,u::SingleFieldFEBasis{<:TrialBasis}
)
  _data = _compute_local_solves(k,u)
  data = lazy_map(transpose,_data)
  return GenericCellField(data,k.trian_out,ReferenceDomain())
end

function Arrays.evaluate!(
  cache,k::LocalOperator,v::CellField
)
  is_test(v) = eltype(get_data(v)) <: AbstractVector{<:Field}
  is_trial(v) = eltype(get_data(v)) <: AbstractMatrix{<:Field}
  u = is_test(v) ? similar_cell_field(v,lazy_map(transpose,get_data(v))) : v
  _data = _compute_local_solves(k,u)
  data = is_trial(v) ? lazy_map(transpose,_data) : _data
  return GenericCellField(data,k.trian_out,ReferenceDomain())
end

function _compute_local_solves(
  k::LocalOperator,u::CellField
)
  cell_coeffs = lazy_map(k.local_map,k.weakform(u))
  if k.collect_coefficients
    cell_coeffs = Arrays.lazy_collect(cell_coeffs)
  end
  v_out = get_fe_basis(k.space_out)
  cell_basis = CellData.get_data(change_domain(v_out,k.trian_out,DomainStyle(v_out)))
  return lazy_map(linear_combination,cell_coeffs,cell_basis)
end

###########################################################################################

struct LocalSolveMap{A} <: Map
  pivot :: A
end

LocalSolveMap() = LocalSolveMap(NoPivot())

function Arrays.evaluate!(cache,k::LocalSolveMap, matvec::Tuple)
  mat, vec = matvec
  evaluate!(cache,k,mat,vec)
end

function Arrays.evaluate!(cache,k::LocalSolveMap, mat, vec)
  f = lu!(mat,k.pivot;check=false)
  @check issuccess(f) "Factorization failed"
  ldiv!(f,vec)
  return vec
end

function Arrays.evaluate!(cache,k::LocalSolveMap, mat::MatrixBlock, vec::MatrixBlock)
  @check size(mat,1) == size(mat,2) == size(vec,1) == 1
  f = lu!(get_array(mat)[1,1],k.pivot;check=false)
  @check issuccess(f) "Factorization failed"
  v = reshape(get_array(vec),size(vec,2))
  for i in eachindex(v)
    ldiv!(f,v[i])
  end
  return v
end

###########################################################################################

struct StaticCondensationMap{A} <: Map
  pivot :: A
end

StaticCondensationMap() = StaticCondensationMap(NoPivot())

function Arrays.evaluate!(cache,k::StaticCondensationMap, matvec::Tuple)
  mat, vec = matvec
  evaluate!(cache,k,mat,vec)
end

function Arrays.evaluate!(cache,k::StaticCondensationMap, mat, vec)
  @check size(mat) == (2,2)
  @check size(vec) == (2,)

  Kii, Kbi, Kib, Kbb = get_array(mat)
  bi, bb = get_array(vec)

  f = lu!(Kii,k.pivot;check=false)
  @check issuccess(f) "Factorization failed"
  ldiv!(f,bi)
  ldiv!(f,Kib)

  mul!(bb,Kbi,bi,-1,1)
  mul!(Kbb,Kbi,Kib,-1,1)

  return Kbb, bb
end

###########################################################################################

struct BackwardStaticCondensationMap{A} <: Map
  pivot :: A
end

BackwardStaticCondensationMap() = BackwardStaticCondensationMap(NoPivot())

function Arrays.evaluate!(cache,k::BackwardStaticCondensationMap, matvec::Tuple, xb)
  mat, vec = matvec
  evaluate!(cache,k,mat,vec,xb)
end

function Arrays.evaluate!(cache,k::BackwardStaticCondensationMap, mat, vec, xb)
  @check size(mat) == (2,2)
  @check size(vec) == (2,)

  Kii, Kbi, Kib, Kbb = get_array(mat)
  bi, bb = get_array(vec)

  f = lu!(Kii, k.pivot; check=false)
  @check issuccess(f) "Factorization failed"

  # Reconstruct interior solution
  mul!(bi, Kib, xb, -1, 1)  # bi = bi - Kib * xb
  ldiv!(f, bi)              # bi = Kii^{-1} * (bi - Kib * xb)

  return bi
end

###########################################################################################
struct HHO_ReconstructionOperatorMap{A} <: Map
  pivot :: A
end

HHO_ReconstructionOperatorMap() = HHO_ReconstructionOperatorMap(NoPivot())

function Arrays.evaluate!(cache::Nothing, k::HHO_ReconstructionOperatorMap, mats::Tuple)
  lhs, rhs = mats
  evaluate!(cache, k, lhs, rhs)
end

function Arrays.evaluate!(cache::Nothing, k::HHO_ReconstructionOperatorMap, lhs, rhs)
  @check size(lhs) == (2,2)
  @check size(rhs) == (2,2)

  App, Aλp, Apλ, _ = get_array(lhs)
  BpΩ, BλΩ, BpΓ, _ = get_array(rhs)

  μT = tr(App)/norm(Apλ)^2
  
  # App = App + μT * Apλ * Aλp
  mul!(App, Apλ, Aλp, μT, 1)
  
  # BpΩ = BpΩ + μT * Apλ * BλΩ
  mul!(BpΩ, Apλ, BλΩ, μT, 1)

  Ainv = lu!(App,k.pivot;check=false)
  @check issuccess(Ainv) "Factorization failed"

  RuΩ = ldiv!(Ainv,BpΩ)
  RuΓ = ldiv!(Ainv,BpΓ)

  return RuΩ, RuΓ
end