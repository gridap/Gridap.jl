
"""
struct ConstantFESpace <: SingleFieldFESpace
  # private fields
end
"""
struct ConstantFESpace{V,T,A,B,C} <: SingleFieldFESpace
model::DiscreteModel
cell_basis::A
cell_dof_basis::B
cell_dof_ids::C
function ConstantFESpace(model;
                       vector_type::Type{V}=Vector{Float64},
                       field_type::Type{T}=Float64) where {V,T}
function setup_cell_reffe(model::DiscreteModel,
  reffe::Tuple{<:ReferenceFEName,Any,Any}; kwargs...)
  basis, reffe_args,reffe_kwargs = reffe
  cell_reffe = ReferenceFE(model,basis,reffe_args...;reffe_kwargs...)
end
reffe=ReferenceFE(lagrangian,T,0)
cell_reffe = setup_cell_reffe(model,reffe)
cell_basis_array=lazy_map(get_shapefuns,cell_reffe)

cell_basis=SingleFieldFEBasis(
  cell_basis_array,
  Triangulation(model),
  TestBasis(),
  ReferenceDomain())

cell_dof_basis_array=lazy_map(get_dof_basis,cell_reffe)
cell_dof_basis=CellDof(cell_dof_basis_array,Triangulation(model),ReferenceDomain())

cell_dof_ids=Fill(Int32(1):Int32(num_components(field_type)),num_cells(model))
A=typeof(cell_basis)
B=typeof(cell_dof_basis)
C=typeof(cell_dof_ids)
new{V,T,A,B,C}(model,
               cell_basis,
               cell_dof_basis,
               cell_dof_ids)
end
end

# Genuine functions
function TrialFESpace(f::ConstantFESpace)
f
end

# Delegated functions
get_triangulation(f::ConstantFESpace) = Triangulation(f.model)

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

function gather_free_and_dirichlet_values!(free_vals,
                                                      dirichlet_vals,
                                                      f::ConstantFESpace,
                                                      cell_vals)
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
