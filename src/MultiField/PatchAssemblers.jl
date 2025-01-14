
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
