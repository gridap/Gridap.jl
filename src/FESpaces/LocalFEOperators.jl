
"""
    struct LocalOperator end

An operator defined by local problems, solved by a local map (e.g. `LocalSolveMap`).

- `LocalOperator(local_map, ptopo, space_to, space_from, lhs, rhs; kwargs...)`: 
  one problem per patch of the `PatchTopology` `ptopo`.
- `LocalOperator(local_map, space_to, lhs, rhs; kwargs...)`: 
  one problem per cell of `trian_out`.

An operator `R` can be applied to

- an `FEFunction` or a `CellField`, returning a `CellField`;
- an FE basis, returning the images of the basis functions;
- an `FESpace` `X`, returning the space `R(X)` that: A) has the dofs of `X` collected onto 
  the local patches, its basis functions are the images of the basis of `X`. 
"""
struct LocalOperator
  local_map    :: Map           # coeffs = evaluate(local_map, lhs, rhs)
  trian_coeffs :: Triangulation # Where the coefficients are defined
  trian_out    :: Triangulation # Where the output field is defined
  space_out    :: FESpace
  weakform     :: Function
  get_out_dof_ids :: Function  # space -> dof ids of `space` on the cells of `trian_out`
  collect_coefficients :: Bool # Collect local coefficients / keep them as lazy maps
end

function LocalOperator(
  local_map::Map,
  ptopo::PatchTopology,
  space_to::FESpace,
  space_from::FESpace,
  lhs::Function,
  rhs::Function;
  space_test::FESpace = space_to,
  space_out::FESpace = space_to,
  trian_out::Triangulation = get_triangulation(space_out),
  collect_coefficients::Bool = true
)
  function weakform(u_from)
    u_to = get_trial_fe_basis(space_to)
    v_to = get_fe_basis(space_test)

    lhs_assem = PatchAssembler(ptopo,space_to,space_test)
    lhs_contr = lhs(u_to,v_to)
    lhs_mats  = assemble_matrix(
      lhs_assem,collect_patch_cell_matrix(lhs_assem,space_to,space_test,lhs_contr)
    )

    rhs_assem = PatchAssembler(ptopo,space_from,space_test)
    rhs_contr = rhs(u_from,v_to)
    if CellData.is_matrix_contribution(rhs_contr)
      rhs_mats = assemble_matrix(
        rhs_assem,collect_patch_cell_matrix(rhs_assem,space_from,space_test,rhs_contr)
      )
    elseif CellData.is_vector_contribution(rhs_contr)
      rhs_mats = assemble_vector(
        rhs_assem,collect_patch_cell_vector(rhs_assem,space_test,rhs_contr)
      )
    else
      @unreachable
    end

    pair_arrays(lhs_mats,rhs_mats)
  end

  trian_coeffs = PatchTriangulation(get_background_model(trian_out),ptopo)

  function get_out_dof_ids(space)
    patch_dof_ids = generate_patch_dof_ids(space,ptopo)
    _move_patch_coeffs(trian_out,trian_coeffs,patch_dof_ids)
  end

  return LocalOperator(
    local_map, trian_coeffs, trian_out, space_out, weakform, get_out_dof_ids, collect_coefficients
  )
end

function LocalOperator(
  local_map::Map,
  space_to::FESpace,
  lhs::Function,
  rhs::Function;
  space_test::FESpace = space_to,
  space_out::FESpace = space_to,
  trian_out::Triangulation = get_triangulation(space_out),
  collect_coefficients::Bool = true
)
  function weakform(u_from)
    u_to = get_trial_fe_basis(space_to)
    v_to = get_fe_basis(space_test)

    lhs_c = lhs(u_to,v_to)
    @check all(t -> t === trian_out, get_domains(lhs_c))
    lhs_mats = get_contribution(lhs_c, trian_out)

    rhs_c = rhs(u_from,v_to)
    @check all(t -> t === trian_out, get_domains(rhs_c))
    rhs_mats = get_contribution(rhs_c, trian_out)

    pair_arrays(lhs_mats,rhs_mats)
  end

  get_out_dof_ids(space) = get_cell_dof_ids(space,trian_out)

  return LocalOperator(
    local_map, trian_out, trian_out, space_out, weakform, get_out_dof_ids, collect_coefficients
  )
end

(P::LocalOperator)(u) = evaluate(P,u)

function Arrays.evaluate!(cache,k::LocalOperator,u)
  _evaluate_local_operator(k,k.space_out,u)
end

function _evaluate_local_operator(k::LocalOperator,space_out::SingleFieldFESpace,u)
  cell_coeffs = _compute_local_solves(k,u)
  return _local_operator_output(k,space_out,u,cell_coeffs,1)
end

function _local_operator_output(
  k::LocalOperator,space_out::SingleFieldFESpace,u::CellField,cell_coeffs,row
)
  _data = _split_local_solves(space_out,k.trian_out,cell_coeffs;row)
  data = is_trial_fe_basis(u) ? lazy_map(transpose,_data) : _data
  domain = DomainStyle(get_fe_basis(space_out))
  return GenericCellField(data,k.trian_out,domain)
end

function _local_operator_output(
  k::LocalOperator,space_out::SingleFieldFESpace,space::SingleFieldFESpace,cell_coeffs,row
)
  cell_shapefuns = _split_local_solves(space_out,k.trian_out,cell_coeffs;row)
  domain = DomainStyle(get_fe_basis(space_out))
  fe_basis = SingleFieldFEBasis(cell_shapefuns,k.trian_out,TestBasis(),domain)
  return PatchFESpace(space,k.trian_out,k.get_out_dof_ids(space),fe_basis)
end

function _compute_local_solves(k::LocalOperator,space::FESpace)
  _compute_local_solves(k,get_trial_fe_basis(space))
end

function _compute_local_solves(k::LocalOperator,u::CellField)
  v = is_test_fe_basis(u) ? as_trial_fe_basis(u) : u
  cell_coeffs = lazy_map(k.local_map,k.weakform(v))
  cell_coeffs = _move_patch_coeffs(k.trian_out,k.trian_coeffs,cell_coeffs)
  if k.collect_coefficients
    cell_coeffs = Arrays.lazy_collect(cell_coeffs)
  end
  return cell_coeffs
end

function _split_local_solves(
  space_out::FESpace, trian_out::Triangulation, cell_coeffs; row::Int = 1
)
  v_out = get_fe_basis(space_out)
  cell_basis = CellData.get_data(change_domain(v_out,trian_out,DomainStyle(v_out)))

  C = testitem(cell_coeffs)
  if C isa AbstractArray
    return lazy_map(linear_combination,cell_coeffs,cell_basis)
  elseif C isa VectorBlock
    nrows = length(C)
    @check row <= nrows
    coeffs = lazy_map(Base.Fix2(Arrays.unsafe_getindex,row), cell_coeffs)
    return lazy_map(linear_combination, coeffs, cell_basis)
  elseif C isa MatrixBlock
    nrows, ncols = size(C)
    @check row <= nrows
    return map(1:ncols) do col
      I = CartesianIndex(row, col)
      coeffs = lazy_map(Base.Fix2(Arrays.unsafe_getindex,I), cell_coeffs)
      lazy_map(linear_combination, coeffs, cell_basis)
    end
  else
    @unreachable
  end
end

function _move_patch_coeffs(trian_out, trian_coeffs, patch_to_coeffs)
  @check trian_out === trian_coeffs
  return patch_to_coeffs
end

function _move_patch_coeffs(trian_out, trian_coeffs::PatchTriangulation, patch_to_coeffs)
  move_contributions(patch_to_coeffs, trian_coeffs, trian_out)
end

function _move_patch_coeffs(
  trian_out::PatchTriangulation, trian_coeffs::PatchTriangulation, patch_to_coeffs
)
  if trian_out === trian_coeffs
    return patch_to_coeffs
  end
  @check trian_out.ptopo === trian_coeffs.ptopo
  tface_to_patch = trian_out.glue.tface_to_patch
  return lazy_map(Reindex(patch_to_coeffs), tface_to_patch)
end
