"""
   abstract type CellDofBasis end

Abstract type that represents a cell array of `Dof`. The main motivation for
its definition is to provide a trait that informs whether the `Dof` entries are
defined for functions in the reference or physical space
"""
abstract type CellDofBasis end

RefStyle(::Type{<:CellDofBasis}) = @abstractmethod
get_array(::CellDofBasis) = @abstractmethod
evaluate(cell_dofs::CellDofBasis,cell_field::CellFieldLike) = @abstractmethod

RefStyle(::T) where T <: CellDofBasis = RefStyle(T)

is_in_ref_space(::Type{T}) where T <:CellDofBasis = get_val_parameter(RefStyle(T))
is_in_ref_space(::T) where T <:CellDofBasis = is_in_ref_space(T)
# is_in_physical_space(a) = !is_in_ref_space(a) # @santiagobadia : not needed, already in CellFieldLike

function test_cell_dof_basis(cf::CellDofBasis,f::CellFieldLike)
  ar = get_array(cf)
  @test isa(ar,AbstractArray)
  a = evaluate(cf,f)
  @test isa(RefStyle(cf),Bool)
end

"""
"""
struct GenericCellDofBasis{R} <: CellDofBasis
   ref_trait::Val{R}
   array::AbstractArray{<:Dof}
end

RefStyle(::Type{<:GenericCellDofBasis{R}}) where R = Val{R}()

get_array(a::GenericCellDofBasis) = a.array

"""
    evaluate(dof_array::CellDofBasis,field_array::AbstractArray)

Evaluates the `CellDofBasis` for the `Field` objects
at the array `field` element by element.

The result is numerically equivalent to

    map(evaluate_dof, dof_array.array, field_array)

but it is described with a more memory-friendly lazy type.
"""
function evaluate(cell_dofs::CellDofBasis,cell_field::CellFieldLike)
 _evaluate_cell_dofs(cell_dofs,cell_field,RefStyle(cell_dofs))
end

function  _evaluate_cell_dofs(cell_dofs,cell_field,ref_trait::Val{true})
  evaluate_dof_array(get_array(cell_dofs),get_array(_to_ref_space(cell_field)),ref_trait)
end

function  _evaluate_cell_dofs(cell_dofs,cell_field,ref_trait::Val{false})
  evaluate_dof_array(get_array(cell_dofs),get_array(_to_physical_space(cell_field)),ref_trait)
end

# @santiagobadia : To be implemented
_to_ref_space(a) = a

_to_physical_space(a) = a
