"""
    abstract type Dof <: Map

Abstract type for a degree of freedom, seen as a linear form over a functional
space (typically a polynomial space). The domain is thus a [`Field`](@ref) set
and the range the scalar set.
"""
abstract type Dof <: Map end

"""
    struct LinearCombinationDofVector{T<:Dof,V,F} <: AbstractVector{T}
      values :: V
      predofs:: F
    end

Type that implements a dof basis (a) as the linear combination of a dof pre-basis
(b). The dofs are first evaluated at the dof pre-basis (b) (field `predofs`) and the
predof values are next mapped to dof basis (a) applying a change of basis (field
`values`).

Fields:

- `values::AbstractMatrix{<:Number}` the matrix of the change from dof basis (b) to (a)
- `predofs::AbstractVector{T}` A type representing dof pre-basis (b), with `T<:Dof`
"""
@ahe struct LinearCombinationDofVector{T,V,F} <: AbstractVector{T}
  values::V
  predofs::F
  function LinearCombinationDofVector(
    values::AbstractMatrix{<:Number},
    predofs::AbstractVector{<:Dof}
  )
    @check size(values,1) == length(predofs) """\n
    Incompatible sizes for performing the linear combination

        linear_combination(values,predofs) = transpose(values)*predofs

    size(values,1) != length(predofs)
    """
    T = eltype(predofs)
    V = typeof(values)
    F = typeof(predofs)
    new{T,V,F}(values,predofs)
  end
end

Base.size(a::LinearCombinationDofVector) = (size(a.values,2),)
Base.IndexStyle(::LinearCombinationDofVector) = IndexLinear()
Base.getindex(::LinearCombinationDofVector{T},::Integer) where T = T()

function linear_combination(a::AbstractMatrix{<:Number}, b::AbstractVector{<:Dof})
  LinearCombinationDofVector(a,b)
end

function return_cache(b::LinearCombinationDofVector,field)
  k = Fields.LinearCombinationMap(:)
  cf = return_cache(b.predofs,field)
  fx = evaluate!(cf,b.predofs,field)
  ck = return_cache(k,fx,transpose(b.values))
  return cf, ck
end

function evaluate!(cache,b::LinearCombinationDofVector,field)
  cf, ck = cache
  k = Fields.LinearCombinationMap(:)
  fx = evaluate!(cf,b.predofs,field)
  return evaluate!(ck,k,fx,transpose(b.values))
end

"""
    struct MappedDofBasis{T<:Dof,MT,BT} <: AbstractVector{T}
      F :: MT
      Σ :: BT
      args
    end

Represents { η = σ∘F : σ ∈ Σ }, evaluated as η(φ) = σ(F(φ,args...)) where

- σ : V -> R  are dofs on V
- F : W  -> V is a map between function spaces

Intended combinations would be:

- Σ ⊂ V* a dof basis in the physical domain and ``F_*`` : V̂ -> V is a pushforward map.
- Σ̂ ⊂ V̂* a dof basis in the reference domain and ``(F_*)``⁻¹ : V -> V̂ is an inverse pushforward map.
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
Base.getindex(::MappedDofBasis{T}, ::Integer) where T = T()

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
    test_dof(dof,field,v; cmp::Function=(==))

Test that the `Dof` interface is properly implemented
for object `dof`. It also checks if the object `dof`
when evaluated at the field `field` returns the same
value as `v`. Comparison is made with the `comp` function.
"""
function test_dof(dof::Dof,field,v;cmp::Function=(==))
  _test_dof(dof,field,v,cmp)
end

"""
    test_dof_array(dofs::AbstractArray{<:Dof}, field, v; cmp::Function=(==))

Like [`test_dof`](@ref) for a DoF basis.
"""
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
