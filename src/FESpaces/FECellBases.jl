
FECellBasisStyle(::Type{T}) where T = Val{false}()

FECellBasisStyle(cell_basis) = FECellBasisStyle(typeof(cell_basis))

FECellBasisStyle(::Type{T}) where T<:CellField = Val(is_basis(T))

"""
"""
function is_a_fe_cell_basis(cell_basis)
  v = FECellBasisStyle(cell_basis)
  get_val_parameter(v)
end

function is_a_fe_cell_basis(::Type)
  @unreachable "is_a_fe_cell_basis cannot be called on types"
end
