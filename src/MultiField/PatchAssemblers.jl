
using Gridap.Geometry: PatchTopology, get_patch_faces
using DataStructures: SortedSet
using BlockArrays

###########################

function FESpaces.get_patch_assembly_ids(space::MultiFieldFESpace,ptopo::Geometry.PatchTopology)
  mfs = MultiFieldStyle(space)
  FESpaces.get_patch_assembly_ids(mfs,space,ptopo)
end

function FESpaces.get_patch_assembly_ids(
  ::ConsecutiveMultiFieldStyle,space::MultiFieldFESpace,ptopo::Geometry.PatchTopology
)
  offsets = compute_field_offsets(space) |> Tuple
  sf_patch_dofs = map(space) do sf
    FESpaces.get_patch_assembly_ids(sf,ptopo)
  end |> Tuple
  mf_patch_dofs = Arrays.append_tables_locally(offsets,sf_patch_dofs)
  return mf_patch_dofs
end

function FESpaces.get_patch_assembly_ids(
  ::BlockMultiFieldStyle{NB,SB,P},space::MultiFieldFESpace,ptopo::Geometry.PatchTopology
) where {NB,SB,P}
  offsets = compute_field_offsets(space) |> Tuple
  sf_patch_dofs = map(space) do sf
    FESpaces.get_patch_assembly_ids(sf,ptopo)
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
  test::MultiFieldFESpace{<:BlockMultiFieldStyle},
)
  NBr, SBr, Pr = get_block_parameters(MultiFieldStyle(test))
  NBc, SBc, Pc = get_block_parameters(MultiFieldStyle(trial))

  block_map = get_block_map(NBr,SBr,Pr,NBc,SBc,Pc)
  block_patch_rows = FESpaces.get_patch_assembly_ids(test,ptopo)
  block_patch_cols = FESpaces.get_patch_assembly_ids(trial,ptopo)
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
  cache,k::FESpaces.LocalOperator,space::MultiFieldFESpace
)
  u = evaluate!(cache,k,get_trial_fe_basis(space))
  v = map(u) do ui
    CellData.similar_cell_field(ui,lazy_map(transpose,CellData.get_data(ui)),get_triangulation(ui),DomainStyle(ui))
  end |> MultiFieldCellField
  return u, v
end

function Arrays.evaluate!(
  cache,k::FESpaces.LocalOperator,u::MultiFieldFEBasisComponent
)
  nfields, fieldid = u.nfields, u.fieldid
  block_fields(fields,::TestBasis) = lazy_map(BlockMap(nfields,fieldid),fields)
  block_fields(fields,::TrialBasis) = lazy_map(BlockMap((1,nfields),fieldid),fields)

  sf = evaluate!(nothing,k,u.single_field)
  data = block_fields(CellData.get_data(sf),BasisStyle(u.single_field))
  return CellData.similar_cell_field(sf,data)
end

@inline function Arrays.evaluate!(
  cache,k::FESpaces.LocalOperator,u::MultiFieldFEFunction
)
  evaluate!(cache,k,u.multi_cell_field)
end

# TODO: The following is a bit of a mess... we would like to return (and allow as inputs)
# a MultiFieldCellField of GenericCellFields. Unfortunately, the changes of domain for 
# blocked cell-data is not implemented. Rather, there are specific functions for the 
# change of domain of MultiFieldFEBasisComponent.
# What tod do? We woudl need to somehow add methods to the `extend`/`pos_neg_data` machinery.
function Arrays.evaluate!(
  cache,k::FESpaces.LocalOperator,u::MultiFieldCellField
)
  block_fields(fields,::TestBasis,i) = lazy_map(BlockMap(nfields,i),fields)
  block_fields(fields,::TrialBasis,i) = lazy_map(BlockMap((1,nfields),i),fields)
  _is_test(v::MultiFieldFEBasisComponent) = isa(BasisStyle(v),TestBasis)
  _is_trial(v::MultiFieldFEBasisComponent) = isa(BasisStyle(v),TrialBasis)
  _is_test(v) = eltype(get_data(v)) <: VectorBlock
  _is_trial(v) = eltype(get_data(v)) <: MatrixBlock

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

  # bstyle = ifelse(is_test, TestBasis(), TrialBasis())
  single_fields = map(1:nfields,mf_data,u) do i, sf_data, u
    sf_data = is_trial ? lazy_map(transpose,sf_data) : sf_data
    # sf_data = is_basis ? block_fields(sf_data,bstyle,i) : sf_data
    # GenericCellField(sf_data,k.trian_out,DomainStyle(u))
    if is_basis
      sf = FESpaces.similar_fe_basis(u.single_field,sf_data,k.trian_out,BasisStyle(u),DomainStyle(u))
      MultiFieldFEBasisComponent(sf,i,nfields)
    else
      GenericCellField(sf_data,k.trian_out,DomainStyle(u))
    end
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

###########################################################################################

struct StaticCondensationOperator <: FEOperator
  full_space :: FESpace
  eliminated_space :: FESpace
  retained_space :: FESpace
  patch_assem :: Assembler
  eliminated_assem :: Assembler
  retained_assem :: Assembler
  full_matvecs
  sc_op
end

function StaticCondensationOperator(
  full_space :: FESpace,
  eliminated_space :: FESpace,
  retained_space :: FESpace,
  patch_assem,
  full_matvecs,
)
  n_elim = num_fields(eliminated_space)
  n_ret = num_fields(retained_space)
  @assert isa(MultiFieldStyle(full_space),BlockMultiFieldStyle{2,(n_elim,n_ret)})

  eliminated_assem = SparseMatrixAssembler(eliminated_space,eliminated_space)
  retained_assem = SparseMatrixAssembler(retained_space,retained_space)

  Asc, bsc = statically_condensed_assembly(retained_assem,patch_assem,full_matvecs)
  sc_op = AffineFEOperator(retained_space,retained_space,Asc,bsc)

  return StaticCondensationOperator(
    full_space,eliminated_space,retained_space,
    patch_assem,eliminated_assem,retained_assem,
    full_matvecs,sc_op
  )
end

function StaticCondensationOperator(
  ptopo,
  full_space :: FESpace,
  eliminated_space :: FESpace,
  retained_space :: FESpace,
  biform :: Function,
  liform :: Function
)
  patch_assem = FESpaces.PatchAssembler(ptopo,full_space,full_space)
  full_matvecs = assemble_matrix_and_vector(biform,liform,patch_assem,full_space,full_space)
  return StaticCondensationOperator(
    full_space,eliminated_space,retained_space,
    patch_assem,full_matvecs
  )
end

function statically_condensed_assembly(retained_assem,patch_assem,full_matvecs)
  sc_matvecs = lazy_map(FESpaces.StaticCondensationMap(),full_matvecs)
  rows = patch_assem.strategy.array.array[2,2].patch_rows
  cols = patch_assem.strategy.array.array[2,2].patch_cols
  data = (([sc_matvecs,],[rows,],[cols,]), ([],[],[]), ([],[]))
  assemble_matrix_and_vector(retained_assem,data)
end

function backward_static_condensation(eliminated_assem,patch_assem,full_matvecs,x_retained)
  rows_elim = patch_assem.strategy.array.array[1,1].patch_rows
  rows_ret = patch_assem.strategy.array.array[2,2].patch_rows

  patch_x_ret = lazy_map(Broadcasting(Reindex(x_retained)),rows_ret)
  patch_x_elim = lazy_map(FESpaces.BackwardStaticCondensationMap(),full_matvecs,patch_x_ret)

  vecdata = ([patch_x_elim,],[rows_elim,])
  assemble_vector(eliminated_assem,vecdata)
end

function backward_static_condensation(op::StaticCondensationOperator,x_retained::AbstractVector)
  backward_static_condensation(op.eliminated_assem,op.patch_assem,op.full_matvecs,x_retained)
end

function backward_static_condensation(op::StaticCondensationOperator,xh_retained)
  x_ret = get_free_dof_values(xh_retained)
  x_elim = backward_static_condensation(op,x_ret)
  return FEFunction(op.eliminated_space,x_elim)
end

FESpaces.get_trial(op::StaticCondensationOperator) = op.full_space
FESpaces.get_test(op::StaticCondensationOperator) = op.full_space
