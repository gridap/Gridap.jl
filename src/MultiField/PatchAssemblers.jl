
using Gridap.Geometry: PatchTopology, get_patch_faces
using DataStructures: SortedSet
using BlockArrays

###########################

function FESpaces.get_patch_assembly_ids(
  space::MultiFieldFESpace, ptopo::Geometry.PatchTopology; kwargs...
)
  mfs = MultiFieldStyle(space)
  FESpaces.get_patch_assembly_ids(mfs,space,ptopo;kwargs...)
end

function FESpaces.get_patch_assembly_ids(
  ::ConsecutiveMultiFieldStyle, space::MultiFieldFESpace, ptopo::Geometry.PatchTopology; kwargs...
)
  offsets = compute_field_offsets(space) |> Tuple
  sf_patch_dofs = map(space) do sf
    FESpaces.get_patch_assembly_ids(sf,ptopo;kwargs...)
  end |> Tuple
  mf_patch_dofs = Arrays.append_tables_locally(offsets,sf_patch_dofs)
  return mf_patch_dofs
end

function FESpaces.get_patch_assembly_ids(
  ::BlockMultiFieldStyle{NB,SB,P}, space::MultiFieldFESpace, ptopo::Geometry.PatchTopology; kwargs...
) where {NB,SB,P}
  offsets = compute_field_offsets(space) |> Tuple
  sf_patch_dofs = map(space) do sf
    FESpaces.get_patch_assembly_ids(sf,ptopo;kwargs...)
  end |> Tuple

  block_ranges = get_block_ranges(NB,SB,P)
  block_patch_dofs = map(block_ranges) do r
    Arrays.append_tables_locally(offsets[r],sf_patch_dofs[r])
  end

  return block_patch_dofs
end

function FESpaces.PatchAssembler(
  ptopo::PatchTopology,
  trial::MultiFieldFESpace{<:BlockMultiFieldStyle},
  test::MultiFieldFESpace{<:BlockMultiFieldStyle};
  kwargs...
)
  NBr, SBr, Pr = get_block_parameters(MultiFieldStyle(test))
  NBc, SBc, Pc = get_block_parameters(MultiFieldStyle(trial))

  block_map = get_block_map(NBr,SBr,Pr,NBc,SBc,Pc)
  block_patch_rows = FESpaces.get_patch_assembly_ids(test,ptopo;kwargs...)
  block_patch_cols = FESpaces.get_patch_assembly_ids(trial,ptopo;kwargs...)
  block_strategies = map(CartesianIndices((NBr,NBc))) do I
    patch_rows = block_patch_rows[I[1]]
    patch_cols = block_patch_cols[I[2]]
    FESpaces.PatchAssemblyStrategy(ptopo,patch_rows,patch_cols)
  end
  strategies = ArrayBlockView(ArrayBlock(block_strategies,fill(true,(NBr,NBc))),block_map)
  rows = map(block_rows -> blockedrange(map(length,block_rows)),zip(block_patch_rows...))
  cols = map(block_cols -> blockedrange(map(length,block_cols)),zip(block_patch_cols...))
  return FESpaces.PatchAssembler(ptopo,strategies,rows,cols)
end

###########################################################################################

function FESpaces.PatchFESpace(space::MultiFieldFESpace,ptopo::PatchTopology)
  spaces = map(space) do sf
    FESpaces.PatchFESpace(sf,ptopo)
  end
  style = MultiFieldStyle(space)
  return MultiFieldFESpace(spaces;style=style)
end

function FESpaces.FESpaceWithoutBCs(space::MultiFieldFESpace)
  spaces = map(FESpaces.FESpaceWithoutBCs,space)
  style = MultiFieldStyle(space)
  return MultiFieldFESpace(spaces;style=style)
end

###########################################################################################

function Arrays.evaluate!(
  cache,k::LocalOperator,space::MultiFieldFESpace
)
  u = evaluate!(cache,k,get_trial_fe_basis(space))
  v = map(u) do ui
    CellData.similar_cell_field(ui,lazy_map(transpose,CellData.get_data(ui)))
  end |> MultiFieldCellField
  return u, v
end

function Arrays.evaluate!(
  cache,k::LocalOperator,u::MultiFieldFEBasisComponent
)
  nfields, fieldid = u.nfields, u.fieldid
  block_fields(fields,::TestBasis) = lazy_map(BlockMap(nfields,fieldid),fields)
  block_fields(fields,::TrialBasis) = lazy_map(BlockMap((1,nfields),fieldid),fields)

  sf = evaluate!(nothing,k,u.single_field)
  return MultiFieldFEBasisComponent(sf,fieldid,nfields)
end

@inline function Arrays.evaluate!(
  cache,k::LocalOperator,u::MultiFieldFEFunction
)
  evaluate!(cache,k,u.multi_cell_field)
end

# TODO: The following is a bit of a mess... we would like to return (and allow as inputs)
# a MultiFieldCellField of GenericCellFields. Unfortunately, the changes of domain for 
# blocked cell-data is not implemented. Rather, there are specific functions for the 
# change of domain of MultiFieldFEBasisComponent.
# What tod do? We woudl need to somehow add methods to the `extend`/`pos_neg_data` machinery.
function Arrays.evaluate!(
  cache,k::LocalOperator,u::MultiFieldCellField
)
  block_fields(fields,::TestBasis,i) = lazy_map(BlockMap(nfields,i),fields)
  block_fields(fields,::TrialBasis,i) = lazy_map(BlockMap((1,nfields),i),fields)
  _is_test(v::MultiFieldFEBasisComponent) = isa(BasisStyle(v),TestBasis)
  _is_trial(v::MultiFieldFEBasisComponent) = isa(BasisStyle(v),TrialBasis)
  _is_test(v) = eltype(get_data(v)) <: VectorBlock
  _is_trial(v) = eltype(get_data(v)) <: MatrixBlock

  domain = DomainStyle(get_fe_basis(k.space_out))
  nfields = num_fields(u)
  is_test  = all(_is_test, u.single_fields)
  is_trial = all(_is_trial, u.single_fields)
  is_basis = is_test || is_trial

  if is_test
    fields = map(1:nfields,u) do i,u
      @assert isa(u,MultiFieldFEBasisComponent)
      sf_data = lazy_map(transpose,CellData.get_data(u.single_field))
      sf = FESpaces.similar_fe_basis(u.single_field,sf_data,get_triangulation(u),TrialBasis(),DomainStyle(u))
      MultiFieldFEBasisComponent(sf,i,nfields)
    end
    v = MultiFieldCellField(fields)
  else
    v = u
  end
  
  mf_data = FESpaces._compute_local_solves(k,v)

  if length(mf_data) != nfields
    @assert isone(length(mf_data))
    return GenericCellField(only(mf_data),k.trian_out,domain)
  end

  single_fields = map(1:nfields,mf_data,u) do i, sf_data, u
    sf_data = is_trial ? lazy_map(transpose,sf_data) : sf_data
    if is_basis
      sf = FESpaces.similar_fe_basis(u.single_field,sf_data,k.trian_out,BasisStyle(u),domain)
      MultiFieldFEBasisComponent(sf,i,nfields)
    else
      GenericCellField(sf_data,k.trian_out,domain)
    end
  end

  return MultiFieldCellField(single_fields)
end
