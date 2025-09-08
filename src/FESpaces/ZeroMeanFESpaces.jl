
"""
    struct ZeroMeanFESpace <: SingleFieldFESpace
      # private fields
    end
"""
struct ZeroMeanFESpace{CA,S} <: SingleFieldFESpace
  space::FESpaceWithConstantFixed{CA,S}
  vol_i::Vector{Float64}
  vol::Float64
end

"""
"""
function ZeroMeanFESpace(space::SingleFieldFESpace,dΩ::Measure)
  _space = FESpaceWithConstantFixed(space,true,num_free_dofs(space))
  vol_i = assemble_vector(v->∫(v)*dΩ,space)
  vol = sum(vol_i)
  ZeroMeanFESpace(_space,vol_i,vol)
end

# Genuine functions

function TrialFESpace(f::ZeroMeanFESpace)
  U = TrialFESpace(f.space)
  ZeroMeanFESpace(U,f.vol_i,f.vol)
end

# function FEFunction(f::ZeroMeanFESpace,free_values)
#   msg = """
#     This function should never be called for ZeroMeanFESpace. Please use 
#     `FEFunction(f::ZeroMeanFESpace,free_values,dirichlet_values)` instead.
#     Reason: 
#       Without the fixed value (i.e the dirichlet_values), we cannot correctly 
#       interpolate the free dofs within the space.
#   """
#   @unreachable msg
# end
# 
function FEFunction(
  f::ZeroMeanFESpace,
  free_values::AbstractVector,
  dirichlet_values::AbstractVector
)
  c = _compute_new_fixedval(
    free_values,
    dirichlet_values,
    f.vol_i,
    f.vol,
    f.space.dof_to_fix
  )
  fv = lazy_map(+,free_values,Fill(c,length(free_values)))
  dv = dirichlet_values .+ c
  FEFunction(f.space,fv,dv)
end

function EvaluationFunction(f::ZeroMeanFESpace,free_values)
  FEFunction(f.space,free_values)
end

function _compute_new_fixedval(fv,dv,vol_i,vol,fixed_dof)
  @assert length(fv) + 1 == length(vol_i)
  @assert length(dv) == 1
  c = zero(eltype(vol_i))
  for i=1:fixed_dof-1
    c += fv[i]*vol_i[i]
  end
  c += first(dv)*vol_i[fixed_dof]
  for i=fixed_dof+1:length(vol_i)
    c += fv[i-1]*vol_i[i]
  end
  c = -c/vol
  return c
end

# This is required, otherwise we end up calling `FEFunction` with a fixed value of zero, 
# which does not properly interpolate the function provided. 
# With this change, we are interpolating in the unconstrained space and then
# substracting the mean.
function interpolate!(object,free_values,fs::ZeroMeanFESpace)
  dirichlet_values = zero_dirichlet_values(fs)
  interpolate_everywhere!(object,free_values,dirichlet_values,fs)
end

# Delegated functions

get_triangulation(f::ZeroMeanFESpace) = get_triangulation(f.space)

ConstraintStyle(::Type{ZeroMeanFESpace{CA,S}}) where {CA,S} = ConstraintStyle(S)

CellField(t::ZeroMeanFESpace,cell_vals) = CellField(t.space,cell_vals)

get_cell_isconstrained(f::ZeroMeanFESpace) = get_cell_isconstrained(f.space)

get_cell_constraints(f::ZeroMeanFESpace) = get_cell_constraints(f.space)

get_dirichlet_dof_values(f::ZeroMeanFESpace) = get_dirichlet_dof_values(f.space)

get_fe_basis(f::ZeroMeanFESpace) = get_fe_basis(f.space)

get_trial_fe_basis(f::ZeroMeanFESpace) = get_trial_fe_basis(f.space)

get_fe_dof_basis(f::ZeroMeanFESpace) = get_fe_dof_basis(f.space)

get_free_dof_ids(f::ZeroMeanFESpace) = get_free_dof_ids(f.space)

zero_free_values(f::ZeroMeanFESpace) = zero_free_values(f.space)

get_vector_type(f::ZeroMeanFESpace) = get_vector_type(f.space)

get_dof_value_type(f::ZeroMeanFESpace) = get_dof_value_type(f.space)

get_cell_dof_ids(f::ZeroMeanFESpace) = get_cell_dof_ids(f.space)

get_dirichlet_dof_ids(f::ZeroMeanFESpace) = get_dirichlet_dof_ids(f.space)

zero_dirichlet_values(f::ZeroMeanFESpace) = zero_dirichlet_values(f.space)

num_dirichlet_tags(f::ZeroMeanFESpace) = num_dirichlet_tags(f.space)

get_dirichlet_dof_tag(f::ZeroMeanFESpace) = get_dirichlet_dof_tag(f.space)

scatter_free_and_dirichlet_values(f::ZeroMeanFESpace,fv,dv) = scatter_free_and_dirichlet_values(f.space,fv,dv)

gather_free_and_dirichlet_values(f::ZeroMeanFESpace,cv) = gather_free_and_dirichlet_values(f.space,cv)

gather_free_and_dirichlet_values!(fv,dv,f::ZeroMeanFESpace,cv) = gather_free_and_dirichlet_values!(fv,dv,f.space,cv)

gather_dirichlet_values(f::ZeroMeanFESpace,cv) = gather_dirichlet_values(f.space,cv)

gather_dirichlet_values!(dv,f::ZeroMeanFESpace,cv) = gather_dirichlet_values!(dv,f.space,cv)

gather_free_values(f::ZeroMeanFESpace,cv) = gather_free_values(f.space,cv)

gather_free_values!(fv,f::ZeroMeanFESpace,cv) = gather_free_values!(fv,f.space,cv)
