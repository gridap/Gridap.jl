
"""
"""
function TestFESpace(args...;kwargs...)
  FESpace(args...;kwargs...)
end

function FESpace(t::Triangulation, reffes; trian=nothing, kwargs...)
  @assert isnothing(trian) || trian === t "The `trian` keyword must be consistant with the given tiangulation $t."
  model = get_active_model(t)
  FESpace(model, reffes; trian=t, kwargs...)
end

function FESpace(
  model::DiscreteModel, 
  reffe::Union{ReferenceFE,Tuple{<:Union{ReferenceFEName,Symbol},Any,Any}}; 
  kwargs...
)
  cell_reffe = ReferenceFE(model,reffe)
  FESpace(model,cell_reffe;kwargs...)
end

"""
    FESpace(model::DiscreteModel, reffe; kwargs...)
    FESpace(trian::Triangulation, reffe; kwargs...)

High level finite element space constructors. The `reffe` can be anything returned by a
[`ReferenceFE`](@ref) constructor, but using one without the `Polytope` argument
is safer, because the polytope(s) is(are) then chosen consistently with `model` or `trian`.

# Keyword Arguments

- `conformity=nothing`: prescribes the space [`Conformity`](@ref), or use `reffe`'s default if `nothing`.

- `dirichlet_tags=Int[]`: tags of entities where Dirichlet boundary conditions are prescribed.

- `dirichlet_masks=nothing`: for Cartesian product finite elements of value type `V::MultiValue`, specifies to which components of `V` do the Dirichlet BC apply, for each tag in `dirichlet_tags`. The type of `dirichlet_masks` is `Nothing` or `Vector{NTuple{N,Bool}}` where [`N = num_indep_components(V)`](@ref num_indep_components).

- `vector_type=nothing`: sets the type of the vector storing the DOF values, e.g. `Vector{Float64}`. Useful to create complex valued function space from a real finite element basis.

- `constraint=nothing`: if set to `:zeromean`, adds a constraint to the space to ensure its functions have zero mean-value.

- `scale_dof=false`: if `true`, rescales the DOF to make them independent from the mesh size. Usefull for mixed-elements using different Piola/physical mappings, that would introduce heterogeneous scaling with meshsize in the systems. By default, the local mesh size is estimated on each physical d-faces/cell using its d-volume.

- `global_meshsize=nothing`: if `scale_dof`, can be used to prescribe a constant homogeneous mesh size estimate, to avoid computing the volume of each physical face owning DOF. This should only be used on very shape-regular and quasi-uniform meshes (e.g. unstretched `CartesianDiscreteModel`).

- `labels=get_face_labeling(model)`: `FaceLabeling` on top of `model` that `dirichlet_tags` refer to.

- `trian=Triangulation(model)`.

# Internal constructors:

    FESpace(model::DiscreteModel, cell_reffe::AbstractArray{<:ReferenceFE}; kwargs...)
    FESpace(model::DiscreteModel, cell_fe::CellFE; kwargs...)
"""
function FESpace(
  model::DiscreteModel,
  cell_reffe::AbstractArray{<:ReferenceFE};
  conformity=nothing,
  trian = Triangulation(model),
  labels = get_face_labeling(model),
  dirichlet_tags=Int[],
  dirichlet_masks=nothing,
  constraint=nothing,
  vector_type=nothing,
  scale_dof=false,
  global_meshsize=nothing
)
  conf = Conformity(testitem(cell_reffe),conformity)
  reg_grid = (num_nodes(model) == num_vertices(model))
  if _use_clagrangian(trian,cell_reffe,conf) && isnothing(constraint) && reg_grid
    return _unsafe_clagrangian(
      cell_reffe,
      Triangulation(model),
      labels,
      vector_type,
      dirichlet_tags,
      dirichlet_masks,
      trian)
  end
  cell_fe = CellFE(model, cell_reffe, conf; scale_dof, global_meshsize)
  return FESpace(
    model,cell_fe; trian, labels, dirichlet_tags, dirichlet_masks, constraint, vector_type
  )
end

function FESpace(
  model::DiscreteModel,
  cell_fe::CellFE;
  trian = Triangulation(model),
  labels = get_face_labeling(model),
  dirichlet_tags=Int[],
  dirichlet_masks=nothing,
  constraint=nothing,
  vector_type=nothing
)
  @assert num_cells(cell_fe) == num_cells(model) """\n
  The number of cells provided in the `cell_fe` argument ($(cell_fe.num_cells) cells)
  does not match the number of cells ($(num_cells(model)) cells) in the provided DiscreteModel.
  """
  conformity = Conformity(cell_fe)
  _vector_type = _get_vector_type(vector_type,cell_fe,trian)
  if isa(conformity,L2Conformity) && isempty(dirichlet_tags)
    F = _DiscontinuousFESpace(_vector_type,trian,cell_fe)
  else
    F = _ConformingFESpace(
      _vector_type, model, labels, cell_fe,
      dirichlet_tags, dirichlet_masks, trian
    )
  end
  V = _add_constraint(F,cell_fe.max_order,constraint)
  return V
end

# Private methods

function _get_vector_type(vector_type,cell_fe,trian)
  if isnothing(vector_type)
    cell_shapefuns, cell_dof_basis = compute_cell_space(cell_fe,trian)
    T = get_dof_value_type(cell_shapefuns,cell_dof_basis)
    return Vector{T}
  else
    @assert vector_type <: AbstractVector """\n
    The (optional) argument vector_type has to be <: AbstractVector.
    """
    return vector_type
  end
end

function _add_constraint(space,order,constraint)
  if isnothing(constraint)
    return space
  elseif constraint == :zeromean
    trian = get_triangulation(space)
    dΩ = Measure(trian,order)
    return ZeroMeanFESpace(space,dΩ)
  end

  @unreachable """\n
    The passed option constraint=$constraint is not valid.
    Valid values for constraint: nothing, :zeromean
  """
end
