
"""
"""
abstract type FESpace <: GridapType end

# Minimal FE interface (used by FEOperator)

"""
"""
function num_free_dofs(f::FESpace)
  @abstractmethod
end

"""
"""
function FEFunction(fe::FESpace, free_values)
  @abstractmethod
end

function EvaluationFunction(fe::FESpace, free_values)
  FEFunction(fe,free_values)
end

"""
"""
function zero_free_values(fs::FESpace)
  @abstractmethod
end

"""
"""
function Base.zero(f::FESpace)
  free_values = zero_free_values(f)
  FEFunction(f,free_values)
end

# Extended FEInterface used by FEOperatorFromTerms and Assemblers

"""
"""
function constraint_style(::Type{<:FESpace})
  @abstractmethod
end

constraint_style(f::T) where T<:FESpace = constraint_style(T)

"""
"""
function has_constraints(::Type{T}) where T <:FESpace
  v = constraint_style(T)
  get_val_parameter(v)
end

has_constraints(f::T) where T<:FESpace = has_constraints(T)

"""
"""
function get_cell_dofs(f::FESpace)
  @abstractmethod
end

function get_cell_dofs(f::FESpace,cellids::AbstractArray)
  reindex(get_cell_dofs(f),cellids)
end

function get_cell_dofs(f::FESpace,cellids::SkeletonPair)
  ids = reindex(get_cell_dofs(f),cellids)
  axs = reindex(get_cell_axes_with_constraints(f),cellids)
  merge_cell_dofs_at_skeleton(ids.left,ids.right,axs.left,axs.right)
end

"""
"""
function get_cell_basis(f::FESpace)
  @abstractmethod
end

"""
"""
function CellData.get_cell_axes(f::FESpace)
  @abstractmethod
end

"""
"""
function get_cell_axes_with_constraints(f::FESpace)
  _get_cell_axes_with_constraints(f,constraint_style(f))
end

function _get_cell_axes_with_constraints(f,::Val{false})
  get_cell_axes(f)
end

function _get_cell_axes_with_constraints(f,::Val{true})
  @abstractmethod
end

function get_cell_constraints(f::FESpace)
  _get_cell_constraints(f,constraint_style(f))
end

function _get_cell_constraints(f,::Val{false})
  identity_constraints(get_cell_axes(f))
end

function _get_cell_constraints(f,::Val{true})
  @abstractmethod
end

function get_cell_constraints(f::FESpace,cellids::AbstractArray)
  reindex(get_cell_constraints(f),cellids)
end

function get_cell_constraints(f::FESpace,cellids::SkeletonPair)
  constr = reindex(get_cell_constraints(f),cellids)
  axsrows = reindex(get_cell_axes_with_constraints(f),cellids)
  axscols = reindex(get_cell_axes(f),cellids)
  merge_cell_constraints_at_skeleton(
    constr.left,constr.right,
    axsrows.left,axsrows.right,
    axscols.left,axscols.right)
end

function get_cell_isconstrained(f::FESpace)
  _get_cell_isconstrained(f,constraint_style(f))
end

function _get_cell_isconstrained(f,::Val{false})
  Fill(false,length(get_cell_dofs(f)))
end

function _get_cell_isconstrained(f,::Val{true})
  @abstractmethod
end

function get_cell_isconstrained(f::FESpace,cellids::AbstractArray)
  reindex(get_cell_isconstrained(f),cellids)
end

function get_cell_isconstrained(f::FESpace,cellids::SkeletonPair)
  isconstr = reindex(get_cell_isconstrained(f),cellids)
  apply((l,r)-> l||r,isconstr.left,isconstr.right)
end

function CellData.attach_constraints_rows(f::FESpace,cellarr,cellids)
  _attach_constraints_rows(f,cellarr,cellids,constraint_style(f))
end

function _attach_constraints_rows(f::FESpace,cellarr,cellids,::Val{false})
  cellarr
end

function _attach_constraints_rows(f::FESpace,cellarr,cellids,::Val{true})
  cellconstr = get_cell_constraints(f,cellids)
  cellmask = get_cell_isconstrained(f,cellids)
  attach_constraints_rows(cellarr,cellconstr,cellmask)
end

function CellData.attach_constraints_cols(f::FESpace,cellarr,cellids)
  _attach_constraints_cols(f,cellarr,cellids,constraint_style(f))
end

function _attach_constraints_cols(f::FESpace,cellarr,cellids,::Val{false})
  cellarr
end

function _attach_constraints_cols(f::FESpace,cellarr,cellids,::Val{true})
  cellconstr = get_cell_constraints(f,cellids)
  cellmask = get_cell_isconstrained(f,cellids)
  attach_constraints_cols(cellarr,cellconstr,cellmask)
end

"""
"""
function test_fe_space(f::FESpace)
  free_values = zero_free_values(f)
  @test length(free_values) == num_free_dofs(f)
  fe_function = FEFunction(f,free_values)
  test_fe_function(fe_function)
  fe_basis = get_cell_basis(f)
  @test isa(has_constraints(f),Bool)
  @test isa(has_constraints(typeof(f)),Bool)
  @test length(get_cell_dofs(f)) == length(fe_basis)
  @test length(get_cell_axes(f)) == length(fe_basis)
  @test length(get_cell_axes_with_constraints(f)) == length(fe_basis)
  @test length(get_cell_constraints(f)) == length(fe_basis)
  @test length(get_cell_isconstrained(f)) == length(fe_basis)
end

function test_fe_space(f::FESpace,matvecdata,matdata,vecdata)
  test_fe_space(f)

  cellmat, cellidsrows, cellidscols = matdata
  cm = attach_constraints_cols(f,cellmat,cellidscols)
  if ! has_constraints(f)
    @test cm === cellmat
  end
  cm = attach_constraints_rows(f,cellmat,cellidsrows)
  if ! has_constraints(f)
    @test cm === cellmat
  end

  cellvec, cellidsrows = vecdata
  cv = attach_constraints_rows(f,cellvec,cellidsrows)
  if ! has_constraints(f)
    @test cv === cellvec
  end

  cellmatvec, cellidsrows, cellidscols = matvecdata
  cmv = attach_constraints_cols(f,cellmatvec,cellidscols)
  if ! has_constraints(f)
    @test cmv === cellmatvec
  end
  cmv = attach_constraints_rows(f,cellmatvec,cellidsrows)
  if ! has_constraints(f)
    @test cmv === cellmatvec
  end

end

