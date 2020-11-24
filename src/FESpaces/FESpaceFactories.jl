
function FESpace(
  model::DiscreteModel,
  cell_fe::CellFE;
  conformity=nothing,
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

  if conformity in (L2Conformity(),:L2) && dirichlet_tags == Int[]
    F = _DiscontinuousFESpace(_vector_type,trian,cell_fe)
  else
    cell_conformity = CellConformity(cell_fe,conformity)
    F = _ConformingFESpace(
      _vector_type,
      model,
      labels,
      cell_fe,
      cell_conformity,
      dirichlet_tags,
      dirichlet_masks)
  end

  V = _add_constraint(F,cell_fe.max_order,constraint)

  V
end

function _add_constraint(F,order,constraint)
  if constraint == nothing
    V = F
  elseif constraint == :zeromean
    trian = get_triangulation(F)
    dΩ = LebesgueMeasure(trian,order)
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
  model::RestrictedDiscreteModel, cell_fe::CellFE;constraint=nothing,kwargs...)
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
  kwargs...)
  trian = Triangulation(model)
  cell_map = get_cell_map(trian)
  if conformity in (nothing, :default, :L2, L2Conformity())
    conf = conformity
  else
    conf = CellConformity(cell_reffe,conformity)
  end
  cell_fe = CellFE(cell_map,cell_reffe)
  FESpace(model,cell_fe;conformity=conf,kwargs...)
end


function FESpace(model::DiscreteModel, reffe::Tuple{Symbol,Any,Any}; kwargs...)
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
