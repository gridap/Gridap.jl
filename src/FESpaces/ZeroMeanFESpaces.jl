
"""
    struct ZeroMeanFESpace <: SingleFieldFESpace
      # private fields
    end
"""
struct ZeroMeanFESpace <: SingleFieldFESpace
  space::FESpaceWithLastDofRemoved
  vol_i::Vector{Float64}
  vol::Float64
end

"""
    ZeroMeanFESpace(
      space::SingleFieldFESpace,
      trian::Triangulation,
      quad::CellQuadrature)
"""
function ZeroMeanFESpace(
  space::SingleFieldFESpace,trian::Triangulation,quad::CellQuadrature)

  _space = FESpaceWithLastDofRemoved(space)
  vol_i, vol = _setup_vols(space,trian,quad)
  ZeroMeanFESpace(_space,vol_i,vol)
end

function _setup_vols(V,trian,quad)
  U = TrialFESpace(V)
  assem = SparseMatrixAssembler(V,U)
  bh = get_cell_basis(V)
  cellvec = integrate(bh,trian,quad)
  cellids = get_cell_id(trian)
  vol_i = assemble_vector(assem,[cellvec],[cellids])
  vol = sum(vol_i)
  (vol_i, vol)
end

# Genuine functions

function TrialFESpace(f::ZeroMeanFESpace)
  U = TrialFESpace(f.space)
  ZeroMeanFESpace(U,f.vol_i,f.vol)
end

function FEFunction(
  f::ZeroMeanFESpace,
  free_values::AbstractVector,
  dirichlet_values::AbstractVector)

  c = _compute_new_fixedval(free_values,dirichlet_values,f.vol_i,f.vol)
  fv = apply(+,free_values,Fill(c,length(free_values)))
  dv = dirichlet_values .+ c
  FEFunction(f.space,fv,dv)
end

function EvaluationFunction(f::ZeroMeanFESpace,free_values)
  FEFunction(f.space,free_values)
end

function _compute_new_fixedval(fv,dv,vol_i,vol)
  @assert length(fv) + 1 == length(vol_i)
  @assert length(dv) == 1
  c = 0.0
  for (i,vi) in enumerate(fv)
    c += vi*vol_i[i]
  end
  c += vol_i[end]*first(dv)
  c = -c/vol
  c
end

# Delegated functions

get_dirichlet_values(f::ZeroMeanFESpace) = get_dirichlet_values(f.space)

get_cell_basis(f::ZeroMeanFESpace) = get_cell_basis(f.space)

get_cell_dof_basis(f::ZeroMeanFESpace) = get_cell_dof_basis(f.space)

num_free_dofs(f::ZeroMeanFESpace) = num_free_dofs(f.space)

zero_free_values(::Type{T},f::ZeroMeanFESpace) where T = zero_free_values(T,f.space)

apply_constraints_matrix_cols(f::ZeroMeanFESpace,cm,cids) = apply_constraints_matrix_cols(f.space,cm,cids)

apply_constraints_matrix_rows(f::ZeroMeanFESpace,cm,cids) = apply_constraints_matrix_rows(f.space,cm,cids)

apply_constraints_vector(f::ZeroMeanFESpace,cm,cids) = apply_constraints_vector(f.space,cm,cids)

get_cell_dofs(f::ZeroMeanFESpace) = get_cell_dofs(f.space)

num_dirichlet_dofs(f::ZeroMeanFESpace) = num_dirichlet_dofs(f.space)

zero_dirichlet_values(f::ZeroMeanFESpace) = zero_dirichlet_values(f.space)

num_dirichlet_tags(f::ZeroMeanFESpace) = num_dirichlet_tags(f.space)

get_dirichlet_dof_tag(f::ZeroMeanFESpace) = get_dirichlet_dof_tag(f.space)

scatter_free_and_dirichlet_values(f::ZeroMeanFESpace,fv,dv) = scatter_free_and_dirichlet_values(f.space,fv,dv)

gather_free_and_dirichlet_values(f::ZeroMeanFESpace,cv) = gather_free_and_dirichlet_values(f.space,cv)

gather_dirichlet_values(f::ZeroMeanFESpace,cv) = gather_dirichlet_values(f.space,cv)

gather_free_values(f::ZeroMeanFESpace,cv) = gather_free_values(f.space,cv)

