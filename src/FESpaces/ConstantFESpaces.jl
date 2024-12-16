
"""
    struct ConstantFESpace <: SingleFieldFESpace
      # private fields
    end

    ConstantFESpace(model::DiscreteModel; vector_type=Vector{Float64}, field_type=Float64)
    ConstantFESpace(trian::Triangulation; vector_type=Vector{Float64}, field_type=Float64)

FESpace that is constant over the provided model/triangulation. Typically used as  
lagrange multipliers. The kwargs `vector_type` and `field_type` are used to specify the
types of the dof-vector and dof-value respectively.
"""
struct ConstantFESpace{V,T,A,B,C} <: SingleFieldFESpace
  trian::Triangulation
  cell_basis::A
  cell_dof_basis::B
  cell_dof_ids::C

  function ConstantFESpace(
    model::DiscreteModel{Dc},
    trian::Triangulation{Dc};
    vector_type::Type{V}=Vector{Float64},
    field_type::Type{T}=Float64
  ) where {Dc,V,T}
    @assert num_cells(model) == num_cells(trian)

    basis, reffe_args, reffe_kwargs = ReferenceFE(lagrangian,T,0)
    cell_reffe = ReferenceFE(model,basis,reffe_args...;reffe_kwargs...)
    cell_basis_array = lazy_map(get_shapefuns,cell_reffe)

    cell_basis = SingleFieldFEBasis(
      cell_basis_array, trian, TestBasis(), ReferenceDomain()
    )
    cell_dof_basis = CellDof(
      lazy_map(get_dof_basis,cell_reffe),trian,ReferenceDomain()
    )
    cell_dof_ids = Fill(Int32(1):Int32(num_indep_components(field_type)),num_cells(trian))

    A = typeof(cell_basis)
    B = typeof(cell_dof_basis)
    C = typeof(cell_dof_ids)
    new{V,T,A,B,C}(trian, cell_basis, cell_dof_basis, cell_dof_ids)
  end
end

function ConstantFESpace(model::DiscreteModel; kwargs...)
  trian = Triangulation(model)
  ConstantFESpace(model,trian; kwargs...)
end

function ConstantFESpace(trian::Triangulation; kwargs...)
  model = get_active_model(trian)
  ConstantFESpace(model,trian; kwargs...)
end

TrialFESpace(f::ConstantFESpace) = f

# Delegated functions
get_triangulation(f::ConstantFESpace) = f.trian

ConstraintStyle(::Type{<:ConstantFESpace}) = UnConstrained()

get_dirichlet_dof_values(f::ConstantFESpace{V}) where V = eltype(V)[]

get_fe_basis(f::ConstantFESpace) = f.cell_basis

get_fe_dof_basis(f::ConstantFESpace) = f.cell_dof_basis

get_free_dof_ids(f::ConstantFESpace) = Base.OneTo(length(f.cell_dof_ids[1]))

get_vector_type(f::ConstantFESpace{V}) where V = V

get_cell_dof_ids(f::ConstantFESpace) = f.cell_dof_ids

get_dirichlet_dof_ids(f::ConstantFESpace) = Base.OneTo(0)

num_dirichlet_tags(f::ConstantFESpace) = 0

get_dirichlet_dof_tag(f::ConstantFESpace) = Int8[]

function scatter_free_and_dirichlet_values(f::ConstantFESpace,fv,dv)
  cell_dof_ids = get_cell_dof_ids(f)
  lazy_map(Broadcasting(PosNegReindex(fv,dv)),cell_dof_ids)
end

function gather_free_and_dirichlet_values!(free_vals, dirichlet_vals, f::ConstantFESpace, cell_vals)
  cell_dofs = get_cell_dof_ids(f)
  cache_vals = array_cache(cell_vals)
  cache_dofs = array_cache(cell_dofs)
  cells = 1:length(cell_vals)

  _free_and_dirichlet_values_fill!(
    free_vals,
    dirichlet_vals,
    cache_vals,
    cache_dofs,
    cell_vals,
    cell_dofs,
    cells)

  (free_vals,dirichlet_vals)
end
