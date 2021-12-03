
function FESpace(t::Triangulation,reffes;trian=nothing, kwargs...)
  @assert trian === nothing
  model = get_active_model(t)
  FESpace(model,reffes;trian=t,kwargs...)
end

function FESpace(
  model::DiscreteModel,
  cell_fe::CellFE;
  trian = Triangulation(model),
  labels = get_face_labeling(model),
  dirichlet_tags=Int[],
  dirichlet_masks=nothing,
  constraint=nothing,
  vector_type=nothing)

  @assert num_cells(cell_fe) == num_cells(model) """\n
  The number of cells provided in the `cell_fe` argument ($(cell_fe.num_cells) cells)
  does not match the number of cells ($(num_cells(model)) cells) in the provided DiscreteModel.
  """
  _vector_type = _get_vector_type(vector_type,cell_fe,trian)
  F = _ConformingFESpace(
      _vector_type,
      model,
      labels,
      cell_fe,
      dirichlet_tags,
      dirichlet_masks,
      trian)
  V = _add_constraint(F,cell_fe.max_order,constraint)
  V
end

function _get_vector_type(vector_type,cell_fe,trian)
  if vector_type == nothing
    cell_shapefuns, cell_dof_basis = compute_cell_space(cell_fe,trian)
    T = get_dof_value_type(cell_shapefuns,cell_dof_basis)
    _vector_type = Vector{T}
  else
    @assert vector_type <: AbstractVector """\n
    The (optional) argument vector_type has to be <: AbstractVector.
    """
    _vector_type = vector_type
  end
  _vector_type
end

function _add_constraint(F,order,constraint)
  if constraint == nothing
    V = F
  elseif constraint == :zeromean
    trian = get_triangulation(F)
    dΩ = Measure(trian,order)
    V = ZeroMeanFESpace(F,dΩ)
  else
    @unreachable """\n
    The passed option constraint=$constraint is not valid.
    Valid values for constraint: nothing, :zeromean
    """
  end
  V
end

function FESpace(
  t::Triangulation,
  cell_reffe::AbstractArray{<:ReferenceFE};
  trian=nothing,
  kwargs...)
  @assert trian === nothing
  # TODO for L2 conformity and no dirichlet conditions
  # no needed to build the active model
  model = get_active_model(t)
  FESpace(model,cell_reffe;trian=t,kwargs...)
end

function FESpace(
  model::DiscreteModel,
  cell_reffe::AbstractArray{<:ReferenceFE};
  conformity=nothing,
  trian = Triangulation(model),
  labels = get_face_labeling(model),
  dirichlet_tags=Int[],
  dirichlet_masks=nothing,
  constraint=nothing,
  vector_type=nothing)

  conf = Conformity(testitem(cell_reffe),conformity)

  if _use_clagrangian(trian,cell_reffe,conf) &&
    constraint === nothing &&
    num_vertices(model) == num_nodes(model)

    V = _unsafe_clagrangian(
      cell_reffe,
      Triangulation(model),
      labels,
      vector_type,
      dirichlet_tags,
      dirichlet_masks,
      trian)
    return V
  end

  cell_fe = CellFE(model,cell_reffe,conf)
  _vector_type = _get_vector_type(vector_type,cell_fe,trian)
  if conformity in (L2Conformity(),:L2) && dirichlet_tags == Int[]
    F = _DiscontinuousFESpace(_vector_type,trian,cell_fe)
    V = _add_constraint(F,cell_fe.max_order,constraint)
  else
    V = FESpace(model,cell_fe;
      trian=trian,
      labels=labels,
      dirichlet_tags=dirichlet_tags,
      dirichlet_masks=dirichlet_masks,
      constraint=constraint,
      vector_type=_vector_type)
  end
  return V
end

function FESpace(model::DiscreteModel,
                 reffe::Tuple{<:ReferenceFEName,Any,Any}; kwargs...)
  basis, reffe_args,reffe_kwargs = reffe
  cell_reffe = ReferenceFE(model,basis,reffe_args...;reffe_kwargs...)
  FESpace(model,cell_reffe;kwargs...)
end

function FESpace(model::DiscreteModel,
                 reffe::ReferenceFE; kwargs...)
  cell_reffe = Fill(reffe,num_cells(model))
  FESpace(model,cell_reffe;kwargs...)
end

"""
"""
function TestFESpace(args...;kwargs...)
  FESpace(args...;kwargs...)
end
