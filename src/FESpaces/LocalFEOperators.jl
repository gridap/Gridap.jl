
"""
    struct LocalOperator end
"""
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

  return LocalOperator(local_map,trian_out,space_out,weakform,collect_coefficients)
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

  return LocalOperator(local_map,trian_out,space_out,weakform,collect_coefficients)
end

(P::LocalOperator)(u) = evaluate(P,u)

function Arrays.evaluate!(
  cache,k::LocalOperator,space::FESpace
)
  u = evaluate!(cache,k,get_trial_fe_basis(space))
  v = CellData.similar_cell_field(u,lazy_map(transpose,CellData.get_data(u)),get_triangulation(u),DomainStyle(u))
  return u, v
end

function Arrays.evaluate!(
  cache,k::LocalOperator,v::SingleFieldFEBasis{<:TestBasis}
)
  u = FESpaces.similar_fe_basis(
    v,lazy_map(transpose,get_data(v)),get_triangulation(v),TrialBasis(),DomainStyle(v)
  )
  data = _compute_local_solves(k,u)
  domain = DomainStyle(get_fe_basis(k.space_out))
  return GenericCellField(data,k.trian_out,domain)
end

function Arrays.evaluate!(
  cache,k::LocalOperator,u::SingleFieldFEBasis{<:TrialBasis}
)
  _data = _compute_local_solves(k,u)
  data = lazy_map(transpose,_data)
  domain = DomainStyle(get_fe_basis(k.space_out))
  return GenericCellField(data,k.trian_out,domain)
end

function Arrays.evaluate!(
  cache,k::LocalOperator,v::CellField
)
  is_test(v) = eltype(get_data(v)) <: AbstractVector{<:Field}
  is_trial(v) = eltype(get_data(v)) <: AbstractMatrix{<:Field}
  u = is_test(v) ? CellData.similar_cell_field(v,lazy_map(transpose,get_data(v))) : v
  _data = _compute_local_solves(k,u)
  data = is_trial(v) ? lazy_map(transpose,_data) : _data
  domain = DomainStyle(get_fe_basis(k.space_out))
  return GenericCellField(data,k.trian_out,domain)
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

  T = eltype(cell_coeffs)
  if T <: AbstractArray
    cell_fields = lazy_map(linear_combination,cell_coeffs,cell_basis)
  elseif T <: Tuple
    nfields = fieldcount(T)
    cell_fields = ntuple(nfields) do i
      coeffs_i = lazy_map(Base.Fix2(getindex,i), cell_coeffs)
      lazy_map(linear_combination, coeffs_i, cell_basis)
    end
  else
    @unreachable
  end

  return cell_fields
end