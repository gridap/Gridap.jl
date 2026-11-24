
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

function FESpaces.PatchFESpace(
  space::MultiFieldFESpace, ptopo::PatchTopology; kwargs...
)
  spaces = map(1:num_fields(space)) do field
    FESpaces.PatchFESpace(space[field],ptopo; kwargs...)
  end
  style = MultiFieldStyle(space)
  return MultiFieldFESpace(spaces;style=style)
end

function FESpaces.PatchFESpace(
  space::MultiFieldFESpace, trian::Triangulation, 
  ptopo::PatchTopology, fe_basis::MultiFieldCellField; kwargs...
)
  @check num_fields(fe_basis) == num_fields(space)
  spaces = map(1:num_fields(space)) do field
    FESpaces.PatchFESpace(space[field],trian,ptopo,fe_basis[field]; kwargs...)
  end
  style = MultiFieldStyle(space)
  return MultiFieldFESpace(spaces;style=style)
end

function FESpaces.PatchFESpace(
  space::MultiFieldFESpace, trian::Triangulation, 
  cell_dof_ids::Vector{<:AbstractArray}, fe_basis::MultiFieldCellField
)
  @check num_fields(fe_basis) == length(cell_dof_ids) == num_fields(space)
  spaces = map(1:num_fields(space)) do field
    FESpaces.PatchFESpace(space[field],trian,cell_dof_ids[field],fe_basis[field])
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
  cache,k::LocalOperator,u::MultiFieldFEBasisComponent
)
  sf = evaluate!(nothing,k,u.single_field)
  return MultiFieldFEBasisComponent(sf,u.fieldid,u.nfields)
end

@inline function Arrays.evaluate!(
  cache,k::LocalOperator,u::MultiFieldFEFunction
)
  evaluate!(cache,k,u.multi_cell_field)
end

# Several outputs: one per field of `space_out`, from a single solve
function FESpaces._evaluate_local_operator(k::LocalOperator,space_out::MultiFieldFESpace,u)
  cell_coeffs = FESpaces._compute_local_solves(k,u)
  return map(1:num_fields(space_out)) do field
    FESpaces._local_operator_output(k,space_out[field],u,cell_coeffs,field)
  end
end

# The output space has one field per field of the input space: the shape functions of field
# `i` are the images of the basis of the i-th input field.
function FESpaces._local_operator_output(
  k::LocalOperator,space_out::SingleFieldFESpace,space::MultiFieldFESpace,cell_coeffs,row
)
  mf_shapefuns = FESpaces._split_local_solves(space_out,k.trian_out,cell_coeffs;row)
  @check length(mf_shapefuns) == num_fields(space)
  domain = DomainStyle(get_fe_basis(space_out))
  fe_basis = MultiFieldCellField(map(mf_shapefuns) do sf_shapefuns
    SingleFieldFEBasis(sf_shapefuns,k.trian_out,TestBasis(),domain)
  end)
  cell_dof_ids = map(k.get_out_dof_ids,space.spaces)
  return FESpaces.PatchFESpace(space,k.trian_out,cell_dof_ids,fe_basis)
end

function FESpaces._local_operator_output(
  k::LocalOperator,space_out::SingleFieldFESpace,u::MultiFieldCellField,cell_coeffs,row
)
  mf_data = FESpaces._split_local_solves(space_out,k.trian_out,cell_coeffs;row)
  domain = DomainStyle(get_fe_basis(space_out))
  is_trial = FESpaces.is_trial_fe_basis(u)
  if !(is_trial || FESpaces.is_test_fe_basis(u))
    return GenericCellField(mf_data,k.trian_out,domain)
  end

  nfields = num_fields(u)
  @check length(mf_data) == nfields
  single_fields = map(1:nfields) do field
    sf_data = is_trial ? lazy_map(transpose,mf_data[field]) : mf_data[field]
    sf = FESpaces.similar_fe_basis(u[field].single_field,sf_data,k.trian_out,BasisStyle(u[field]),domain)
    MultiFieldFEBasisComponent(sf,field,nfields)
  end
  return MultiFieldCellField(single_fields)
end
