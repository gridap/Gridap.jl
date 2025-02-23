
using Gridap.Geometry: PatchTopology, get_patch_faces
using DataStructures: SortedSet
using BlockArrays

struct BlockAssemblyStrategyMap{S,A}
  blocks::A
  function BlockAssemblyStrategyMap{S}(strategies) where S
    blocks = map(AssemblyStrategyMap{S},strategies)
    A = typeof(blocks)
    new{S,A}(blocks)
  end
  function BlockAssemblyStrategyMap{S}(strategies,block_map) where S
    blocks = map(block_map) do I
      FESpaces.AssemblyStrategyMap{S}(strategies[I])
    end
    A = typeof(blocks)
    new{S,A}(blocks)
  end
end

function FESpaces.AssemblyStrategyMap{S}(strategies::ArrayBlockView) where S
  BlockAssemblyStrategyMap{S}(strategies.array.array,strategies.block_map)
end

function Arrays.return_cache(k::BlockAssemblyStrategyMap,ids::ArrayBlock,patch)
  si = testitem(k.blocks)
  fi = testitem(ids)
  ci = return_cache(si,fi,patch)
  gi = evaluate!(ci,si,fi,patch)
  b = Array{typeof(ci),ndims(ids)}(undef,size(ids))
  for i in eachindex(ids.array)
    if ids.touched[i]
      ki = return_cache(si,ids.array[i])
      b[i] = return_cache(si,ids.array[i])
    end
  end
  array = Array{typeof(gi),ndims(ids)}(undef,size(ids))
  ArrayBlock(array,ids.touched), b
end

function Arrays.evaluate!(cache,k::BlockAssemblyStrategyMap{:rows},ids::ArrayBlock,patch)
  a,b = cache
  for i in eachindex(ids.array)
    if ids.touched[i]
      a.array[i] = evaluate!(b[i],k.blocks[i,1],ids.array[i],patch)
    end
  end
  a
end

function Arrays.evaluate!(cache,k::BlockAssemblyStrategyMap{:cols},ids::ArrayBlock,patch)
  a,b = cache
  for i in eachindex(ids.array)
    if ids.touched[i]
      a.array[i] = evaluate!(b[i],k.blocks[1,i],ids.array[i],patch)
    end
  end
  a
end

###########################

function FESpaces.get_patch_dofs(space::MultiFieldFESpace,ptopo::Geometry.PatchTopology)
  mfs = MultiFieldStyle(space)
  FESpaces.get_patch_dofs(mfs,space,ptopo)
end

function FESpaces.get_patch_dofs(
  ::ConsecutiveMultiFieldStyle,space::MultiFieldFESpace,ptopo::Geometry.PatchTopology
)
  offsets = compute_field_offsets(space) |> Tuple
  sf_patch_dofs = map(space) do sf
    FESpaces.get_patch_dofs(sf,ptopo)
  end |> Tuple
  mf_patch_dofs = Arrays.append_tables_locally(offsets,sf_patch_dofs)
  return mf_patch_dofs
end

function FESpaces.get_patch_dofs(
  ::BlockMultiFieldStyle{NB,SB,P},space::MultiFieldFESpace,ptopo::Geometry.PatchTopology
) where {NB,SB,P}
  offsets = compute_field_offsets(space) |> Tuple
  sf_patch_dofs = map(space) do sf
    FESpaces.get_patch_dofs(sf,ptopo)
  end |> Tuple

  block_ranges = get_block_ranges(NB,SB,P)
  block_patch_dofs = map(block_ranges) do r
    Arrays.append_tables_locally(offsets[r],sf_patch_dofs[r])
  end

  return block_patch_dofs
end

function FESpaces.PatchAssembler(
  ptopo::PatchTopology,trial::MultiFieldFESpace{MS},test::MultiFieldFESpace{MS}
) where MS <: BlockMultiFieldStyle{NB,SB,P} where {NB,SB,P}
  NV = length(test.spaces)
  block_map = get_block_map(NB,NV,SB,P)
  block_patch_rows = FESpaces.get_patch_dofs(test,ptopo)
  block_patch_cols = FESpaces.get_patch_dofs(trial,ptopo)
  block_strategies = map(CartesianIndices((NB,NB))) do I
    patch_rows = block_patch_rows[I[1]]
    patch_cols = block_patch_cols[I[2]]
    FESpaces.PatchAssemblyStrategy(ptopo,patch_rows,patch_cols)
  end
  strategies = ArrayBlockView(ArrayBlock(block_strategies,fill(true,(NB,NB))),block_map)
  rows = map(block_rows -> blockedrange(map(length,block_rows)),zip(block_patch_rows...))
  cols = map(block_cols -> blockedrange(map(length,block_cols)),zip(block_patch_cols...))
  return FESpaces.PatchAssembler(ptopo,strategies,rows,cols)
end

###########################################################################################

function Arrays.evaluate!(
  cache,k::FESpaces.LocalOperator,u::MultiFieldFEBasisComponent
)
  nfields, fieldid = u.nfields, u.fieldid
  block_fields(fields,::TestBasis) = lazy_map(BlockMap(nfields,fieldid),fields)
  block_fields(fields,::TrialBasis) = lazy_map(BlockMap((1,nfields),fieldid),fields)

  sf = evaluate!(nothing,k,u.single_field)
  data = block_fields(CellData.get_data(sf),BasisStyle(u.single_field))
  return CellData.similar_cell_field(sf,data,get_triangulation(sf),DomainStyle(sf))
end

function Arrays.evaluate!(
  cache,k::FESpaces.LocalOperator,u::MultiFieldFEFunction
)
  evaluate!(cache,k,u.multi_cell_field)
end

function Arrays.evaluate!(
  cache,k::FESpaces.LocalOperator,u::MultiFieldCellField
)
  block_fields(fields,::TestBasis,i) = lazy_map(BlockMap(nfields,i),fields)
  block_fields(fields,::TrialBasis,i) = lazy_map(BlockMap((1,nfields),i),fields)

  is_basis = all(f -> isa(f,MultiFieldFEBasisComponent), u.single_fields)
  is_test  = is_basis && all(f -> isa(BasisStyle(f),TestBasis), u.single_fields)
  is_trial = is_basis && !is_test
  bstyle = ifelse(is_test, TestBasis(), TrialBasis())
  v = ifelse(is_test, MultiFieldCellField(map(transpose,u)), u)
  
  mf_data = FESpaces._compute_local_solves(k,v)

  nfields = num_fields(u)
  single_fields = map(1:nfields,mf_data,u) do i, sf_data, sf
    sf_data = ifelse(is_trial, lazy_map(transpose,sf_data), sf_data)
    sf_data = ifelse(is_basis, block_fields(sf_data,bstyle,i), sf_data)
    GenericCellField(sf_data,k.trian_out,DomainStyle(sf))
  end

  return MultiFieldCellField(single_fields)
end

function FESpaces._compute_local_solves(
  k::FESpaces.LocalOperator,u::MultiFieldCellField
)
  nfields = num_fields(u)
  cell_coeffs = lazy_map(k.local_map,k.weakform(u))
  if k.collect_coefficients
    cell_coeffs = Arrays.lazy_collect(cell_coeffs)
  end
  v_out = get_fe_basis(k.space_out)
  cell_basis = CellData.get_data(change_domain(v_out,k.trian_out,DomainStyle(v_out)))
  cell_fields = map(1:nfields) do i
    coeffs_i = map(x -> getindex(x,i), cell_coeffs)
    lazy_map(linear_combination,coeffs_i,cell_basis)
  end
  return cell_fields
end
