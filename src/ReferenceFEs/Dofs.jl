abstract type Dof <: Map end

"""
    struct LinearCombinationDofVector{T} <: AbstractVector{Dof}
      values::Matrix{T}
      dofs::AbstractVector{<:Dof}
    end

Type that implements a dof basis (a) as the linear combination of a dof basis
(b). The dofs are first evaluated at dof basis (b) (field `dofs`) and the
dof values are next mapped to dof basis (a) applying a change of basis (field
`values`).

Fields:

- `values::Matrix{T}` the matrix of the change from dof basis (b) to (a)
- `dofs::AbstractVector{<:Dof}` A type representing dof basis (b)
"""
struct LinearCombinationDofVector{T,V,F} <: AbstractVector{T}
  values::V
  dofs::F
  function LinearCombinationDofVector(
    values::AbstractMatrix{<:Number},
    dofs::AbstractVector{<:Dof}
  )
    T = eltype(dofs)
    V = typeof(values)
    F = typeof(dofs)
    new{T,V,F}(values,dofs)
  end
end

@inline Base.size(a::LinearCombinationDofVector) = size(a.values,2)
@inline Base.IndexStyle(::LinearCombinationDofVector) = IndexLinear()
@inline Base.getindex(::LinearCombinationDofVector{T},::Integer) where T = T()

function linear_combination(a::AbstractMatrix{<:Number}, b::AbstractVector{<:Dof})
  LinearCombinationDofVector(a,b)
end

function return_cache(b::LinearCombinationDofVector,field)
  k = LinearCombinationMap(:)
  cf = return_cache(b.dofs,field)
  fx = evaluate!(cf,b.dofs,field)
  ck = return_cache(k,b.values,fx)
  return cf, ck
end

function evaluate!(cache,b::LinearCombinationDofVector,field)
  cf, ck = cache
  k = LinearCombinationMap(:)
  fx = evaluate!(cf,b.dofs,field)
  return evaluate!(ck,k,b.values,fx)
end

"""
    test_dof(dof,field,v;cmp::Function=(==))

Test that the `Dof` interface is properly implemented
for object `dof`. It also checks if the object `dof`
when evaluated at the field `field` returns the same
value as `v`. Comparison is made with the `comp` function.
"""
function test_dof(dof::Dof,field,v;cmp::Function=(==))
  _test_dof(dof,field,v,cmp)
end

function test_dof_array(dof::AbstractArray{<:Dof},field,v;cmp::Function=(==))
  _test_dof(dof,field,v,cmp)
end

function _test_dof(dof,field,v,cmp)
  if isa(dof,Dof)
    test_map(v,dof,field;cmp=cmp)
  end
  r = evaluate(dof,field)
  @test cmp(r,v)
  @test typeof(r) == return_type(dof,field)
end
