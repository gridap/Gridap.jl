"""
   abstract type CellDofBasis end

Abstract type that represents a cell array of `Dof`. The main motivation for
its definition is to provide a trait that informs whether the `Dof` entries are
defined for functions in the reference or physical space
"""
abstract type CellDofBasis end

struct GenericCellDofBasis{R} <: CellDofBasis
   ref_trait::Val{R}
   array::AbstractArray{<:Dof}

   GenericCellDofBasis(R,array) = new{R}(Val{R}(),array)

end

RefTrait(::Type{<:GenericCellDofBasis{R}}) where R = Val{R}()
RefTrait(a::GenericCellDofBasis) = RefTrait(typeof(a))

get_array(a::GenericCellDofBasis) = a.array

"""
    evaluate(dof_array::CellDofBasis,field_array::AbstractArray)

Evaluates the `CellDofBasis` for the `Field` objects
at the array `field` element by element.

The result is numerically equivalent to

    map(evaluate_dof, dof_array.array, field_array)

but it is described with a more memory-friendly lazy type.
"""
function evaluate(cell_dofs::CellDofBasis,cell_field) #::CellFieldLike)
 _evaluate_cell_dofs(cell_dofs,cell_field,RefTrait(cell_dofs))
end

function  _evaluate_cell_dofs(cell_dofs,cell_field,ref_trait::Val{true})
  evaluate_dof_array(get_array(cell_dofs),get_array(_to_ref_space(cell_field)),ref_trait)
end

function  _evaluate_cell_dofs(cell_dofs,cell_field,ref_trait::Val{false})
  evaluate_dof_array(get_array(cell_dofs),get_array(_to_physical_space(cell_field)),ref_trait)
end

# function evaluate_dof_array(dof_array::AbstractArray,field_array::AbstractArray,::Val{true})
#   k = DofEval()
#   apply(k,dof_array,field)
# end
#
# function evaluate_dof_array(dof_array::AbstractArray,field_array::AbstractArray,::Val{false})
#   apply(dof_array,field_array)
# end

_to_ref_space(a) = a

_to_physical_space(a) = a
