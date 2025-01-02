abstract type Dof <: Map end

"""
    struct LinearCombinationDofVector{T<:Dof,V,F} <: AbstractVector{T}
      values :: V
      dofs   :: F
    end

Type that implements a dof basis (a) as the linear combination of a dof basis
(b). The dofs are first evaluated at dof basis (b) (field `dofs`) and the
dof values are next mapped to dof basis (a) applying a change of basis (field
`values`).

Fields:

- `values::AbstractMatrix{<:Number}` the matrix of the change from dof basis (b) to (a)
- `dofs::AbstractVector{T}` A type representing dof basis (b), with `T<:Dof`
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

Base.size(a::LinearCombinationDofVector) = size(a.values,2)
Base.IndexStyle(::LinearCombinationDofVector) = IndexLinear()
Base.getindex(::LinearCombinationDofVector{T},::Integer) where T = T()

function linear_combination(a::AbstractMatrix{<:Number}, b::AbstractVector{<:Dof})
  LinearCombinationDofVector(a,b)
end

function return_cache(b::LinearCombinationDofVector,field)
  k = Fields.LinearCombinationMap(:)
  cf = return_cache(b.dofs,field)
  fx = evaluate!(cf,b.dofs,field)
  ck = return_cache(k,b.values,fx)
  return cf, ck
end

function evaluate!(cache,b::LinearCombinationDofVector,field)
  cf, ck = cache
  k = Fields.LinearCombinationMap(:)
  fx = evaluate!(cf,b.dofs,field)
  return evaluate!(ck,k,b.values,fx)
end

"""
    struct MappedDofBasis{T<:Dof,MT,BT} <: AbstractVector{T}
      F :: MT
      σ :: BT
      args
    end

Represents η = σ∘F, evaluated as η(φ) = σ(F(φ,args...))

  - σ : V* -> R is a dof basis
  - F : W  -> V is a map between function spaces

Intended combinations would be: 

- σ : V* -> R dof basis in the physical domain and F* : V̂ -> V is a pushforward map.
- ̂σ : V̂* -> R dof basis in the reference domain and (F*)^-1 : V -> V̂ is an inverse pushforward map.

"""
struct MappedDofBasis{T<:Dof,MT,BT,A} <: AbstractVector{T}
  F    :: MT
  dofs :: BT
  args :: A

  function MappedDofBasis(F::Map, dofs::AbstractVector{<:Dof}, args...)
    T  = eltype(dofs)
    MT = typeof(F)
    BT = typeof(dofs)
    A  = typeof(args)
    new{T,MT,BT,A}(F,dofs,args)
  end
end

Base.size(b::MappedDofBasis) = size(b.dofs)
Base.IndexStyle(::MappedDofBasis) = IndexLinear()
Base.getindex(b::MappedDofBasis, i) = T()

function Arrays.return_cache(b::MappedDofBasis, fields)
  f_cache = return_cache(b.F,fields,b.args...)
  ffields = evaluate!(f_cache,b.F,fields,b.args...)
  dofs_cache = return_cache(b.dofs,ffields)
  return f_cache, dofs_cache
end

function Arrays.evaluate!(cache, b::MappedDofBasis, fields)
  f_cache, dofs_cache = cache
  ffields = evaluate!(f_cache,b.F,fields,b.args...)
  evaluate!(dofs_cache,b.dofs,ffields)
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
