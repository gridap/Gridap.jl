
function FESpace(
  model::DiscreteModel,
  cell_reffes::AbstractArray{<:ReferenceFE};
  labels = get_face_labeling(model),
  conformity=Conformity(first(cell_reffes)),
  dirichlet_tags=Int[],
  dirichlet_masks=nothing,
  constraint=nothing,
  vectortype::Union{Nothing,Type}=nothing)

  @assert length(cell_reffes) == num_cells(model) """\n
  The length of the vector provided in the `cell_reffes` argument ($(length(cell_reffes)) entries)
  does not match the number of cells ($(num_cells(model)) cells) in the provided DiscreteModel.
  """

  trian = get_triangulation(model)
  shapefuns, dof_basis = compute_cell_space(cell_reffes,trian)

  if vectortype == nothing
    T = get_dof_value_type(shapefuns,dof_basis)
    _vector_type = Vector{T}
  else
    @assert vectortype <: AbstractVector """\n
    """
    _vector_type = vectortype
  end

  _conformity = Conformity(first(cell_reffes),conformity)

  if _conformity == L2Conformity() && dirichlet_tags == Int[]

    F = _DiscontinuousFESpace(
      vector_type,
      trian,
      cell_reffes,
      cell_shapefuns,
      cell_dof_basis)

  else

    F = _ConformingFESpace(
      _vector_type,
      model,
      labels,
      cell_reffes,
      shapefuns,
      dof_basis,
      _conformity,
      dirichlet_tags,
      dirichlet_masks)
  end


  if constraint == nothing
    V = F
  elseif constraint == :zeromean
    @notimplemented "zeromean option not yet implemented"
  else
    @unreachable """\n
    The passed option constraint=$constraint is not valid.
    Valid values for constraint: nothing, :zeromean
    """
  end

  V
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
