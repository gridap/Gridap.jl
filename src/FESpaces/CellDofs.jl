"""
   abstract type CellDofBasis end

Abstract type that represents a cell array of `Dof`. The main motivation for
its definition is to provide a trait that informs whether the `Dof` entries are
defined for functions in the reference or physical space
"""
abstract type CellDofBasis end

# RefTrait{::Type{<:CellDofBasis}} = Val{true}()

struct GenericCellDofBasis{R} <: CellDofBasis
   ref_trait::Val{R}
   array::AbstractArray{<:Dof}
end

RefTrait(::Type{<:GenericCellDofBasis{R}}) where R = Val{R}()
