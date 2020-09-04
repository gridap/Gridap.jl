"""
   abstract type CellDof end

Abstract type that represents a cell array of `Dof`. The main motivation for
its definition is to provide a trait that informs whether the `Dof` entries are
defined for functions in the reference or physical space
"""
abstract type CellDof end

get_array(::CellDof) = @abstractmethod

function test_cell_dof_basis(cf::CellDof,f::CellField)
  a = evaluate(cf,f)
  _ = collect(a)
end

"""
    evaluate(dof_array::CellDof,field_array::AbstractArray)

Evaluates the `CellDof` for the `Field` objects
at the array `field` element by element.

The result is numerically equivalent to

    map(evaluate_dof, dof_array.array, field_array)

but it is described with a more memory-friendly lazy type.
"""
function evaluate(cell_dofs::CellDof,cell_field::CellField)
  ReferenceFEs.evaluate_dof_array(get_array(cell_dofs),get_array(cell_field))
end

"""
Functor-like evaluation for CellDof
"""
(a::CellDof)(f::CellField) = evaluate(a,f)

"""
"""
struct GenericCellDof <: CellDof
   array::AbstractArray{<:Dof}
end

get_array(a::GenericCellDof) = a.array

function reindex(cf::CellDof,a::AbstractVector)
  GenericCellDof(reindex(get_array(cf),a))
end

"""
    struct MappedCellDof
      σ_ref::CellDof
      ϕ::CellField
    end

Basis `σ`  defined as `σ(f) = σ_ref(f∘ϕ) `
"""
struct MappedCellDof
  σ_ref::CellDof
  ϕ::CellField
end

function evaluate(a::MappedCellDof,f::CellField)
  a.σ_ref(f∘a.ϕ)
end

