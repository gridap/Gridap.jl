
function FESpace(
  model::DiscreteModel,
  cell_fe::CellFE;
  labels = get_face_labeling(model),
  dirichlet_tags=Int[],
  dirichlet_masks=nothing,
  constraint=nothing,
  vector_type::Union{Nothing,Type}=nothing)

  @assert num_cells(cell_fe) == num_cells(model) """\n
  The number of cells provided in the `cell_fe` argument ($(cell_fe.num_cells) cells)
  does not match the number of cells ($(num_cells(model)) cells) in the provided DiscreteModel.
  """
  trian = Triangulation(model)
  _vector_type = _get_vector_type(vector_type,cell_fe,trian)
  F = _ConformingFESpace(
      _vector_type,
      model,
      labels,
      cell_fe,
      dirichlet_tags,
      dirichlet_masks)
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
  model::RestrictedDiscreteModel,
  cell_reffe::AbstractArray{<:ReferenceFE};
  conformity=nothing,
  constraint=nothing,kwargs...)
  model_portion = model.model
  conf = Conformity(testitem(cell_reffe),conformity)
  cell_fe = CellFE(model,cell_reffe,conf)
  V_portion = FESpace(model_portion,cell_fe;constraint=nothing,kwargs...)
  F = ExtendedFESpace(V_portion,model)
  V = _add_constraint(F,cell_fe.max_order,constraint)
  V
end

function FESpace(
  model::RestrictedDiscreteModel,
  cell_fe::CellFE;
  constraint=nothing,kwargs...)
  model_portion = model.model
  V_portion = FESpace(model_portion,cell_fe;constraint=nothing,kwargs...)
  F = ExtendedFESpace(V_portion,model)
  V = _add_constraint(F,cell_fe.max_order,constraint)
  V
end

function FESpace(
  model::DiscreteModel,
  cell_reffe::AbstractArray{<:ReferenceFE};
  conformity=nothing,
  labels = get_face_labeling(model),
  dirichlet_tags=Int[],
  dirichlet_masks=nothing,
  constraint=nothing,
  vector_type::Union{Nothing,Type}=nothing)

  trian = Triangulation(model)
  conf = Conformity(testitem(cell_reffe),conformity)

  if _use_clagrangian(trian,cell_reffe,conf) && constraint === nothing
    V = _unsafe_clagrangian(
      cell_reffe,
      trian,
      labels,
      vector_type,
      dirichlet_tags,
      dirichlet_masks)
    return V
  end

  cell_fe = CellFE(model,cell_reffe,conf)
  _vector_type = _get_vector_type(vector_type,cell_fe,trian)
  if conformity in (L2Conformity(),:L2) && dirichlet_tags == Int[]
    F = _DiscontinuousFESpace(_vector_type,trian,cell_fe)
    V = _add_constraint(F,cell_fe.max_order,constraint)
  else
    V = FESpace(model,cell_fe;
      labels=labels,
      dirichlet_tags=dirichlet_tags,
      dirichlet_masks=dirichlet_masks,
      constraint=constraint,
      vector_type=_vector_type)
  end
  return V
end


function FESpace(model::DiscreteModel, reffe::Tuple{<:ReferenceFEName,Any,Any}; kwargs...)
  basis, reffe_args,reffe_kwargs = reffe
  cell_reffe = ReferenceFE(model,basis,reffe_args...;reffe_kwargs...)
  FESpace(model,cell_reffe;kwargs...)
end

function FESpace(model::DiscreteModel, reffe::ReferenceFE; kwargs...)
  cell_reffe = Fill(reffe,num_cells(model))
  FESpace(model,cell_reffe;kwargs...)
end

"""
"""
function TestFESpace(args...;kwargs...)
  FESpace(args...;kwargs...)
end
