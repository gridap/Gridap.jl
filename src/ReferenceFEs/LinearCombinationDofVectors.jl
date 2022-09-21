"""
    struct LinearCombinationDofVector{T} <: AbstractVector{Dof}
      change_of_basis::Matrix{T}
      dof_basis::AbstractVector{<:Dof}
    end

Type that implements a dof basis (a) as the linear combination of a dof basis
(b). The dofs are first evaluated at dof basis (b) (field `dof_basis`) and the
dof values are next mapped to dof basis (a) applying a change of basis (field
`change_of_basis`).

Fields:

- `change_of_basis::Matrix{T}` the matrix of the change from dof basis (b) to (a)
- `dof_basis::AbstractVector{<:Dof}` A type representing dof basis (b)
"""
struct LinearCombinationDofVector{T} <: AbstractVector{Dof}
  change_of_basis::Matrix{T}
  dof_basis::AbstractVector{<:Dof}
end

@inline Base.size(a::LinearCombinationDofVector) = size(a.dof_basis)
@inline Base.axes(a::LinearCombinationDofVector) = axes(a.dof_basis)
@inline Base.getindex(a::LinearCombinationDofVector,i::Integer) = getindex(a.dof_basis,i)
@inline Base.IndexStyle(::LinearCombinationDofVector) = IndexLinear()

function linear_combination(a::AbstractMatrix{<:Number},
                            b::AbstractVector{<:Dof})
  LinearCombinationDofVector(a,b)
end

function linear_combination(a::LinearCombinationDofVector{T},
                            b::AbstractVector{<:Dof}) where T
  linear_combination(a.change_of_basis,b)
end

function return_cache(b::LinearCombinationDofVector,field)
  c, cf = return_cache(b.dof_basis,field)
  c, cf, return_cache(*,b.change_of_basis,c)
end

@inline function evaluate!(cache,b::LinearCombinationDofVector,field)
  c, cf, cc = cache
  vals = evaluate!(cache,b.dof_basis,field)
  evaluate!(cc,*,b.change_of_basis,vals)
end