
function FESpace(
  model::DiscreteModel,
  reffes::AbstractArray{<:ReferenceFE};
  labels = get_face_labeling(model),
  dof_space=:reference,
  conformity::Conformity=get_default_conformity(first(reffes)),
  dirichlet_tags=Int[],
  dirichlet_masks=nothing,
  constraint=nothing,
  vectortype::Union{Nothing,Type}=nothing)

  @assert length(reffes) == num_cells(model) """\n
  The length of the vector provided in the `reffes` argument ($(length(reffes)) entries)
  does not match the number of cells ($(num_cells(model)) cells) in the provided DiscreteModel.
  """

  if dof_space == :reference
    domain_style = ReferenceDomain()
  elseif dof_space == :physical
    domain_style = PhysicalDomain()
  else
    @unreachable """\n
    The passed option dof_space=$dof_space is not valid.
    Valid values for dof_space: :reference, :physical
    """
  end

  trian = get_triangulation(model)
  shapefuns, dof_basis = compute_cell_space(reffes,trian,domain_style)

  if vectortype == nothing
    T = get_dof_value_type(shapefuns,dof_basis)
    _vector_type = Vector{T}
  else
    @assert vectortype <: AbstractVector """\n
    """
    _vector_type = vectortype
  end

  F = _ConformingFESpace(
    _vector_type,
    model,
    labels,
    reffes,
    shapefuns,
    dof_basis,
    conformity,
    dirichlet_tags,
    dirichlet_masks)

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

function FESpace(model::DiscreteModel, reffe::Tuple{Symbol,Any}; kwargs...)
  basis, reffe_kwargs = reffe
  cell_reffe = ReferenceFE(model,basis;reffe_kwargs...)
  FESpace(model,cell_reffe;kwargs...)
end

"""
"""
function TestFESpace(args...;kwargs...)
  FESpace(args...;kwargs...)
end
