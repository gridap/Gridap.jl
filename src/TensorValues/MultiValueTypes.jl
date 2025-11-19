###############################################################
# MultiValue Type
###############################################################

"""
    MultiValue{S,T,N,L} <: Number

Abstract type representing a multi-dimensional number value. The parameters are analog to that of StaticArrays.jl:
- `S` is a Tuple type holding the size of the tensor, e.g. Tuple{3} for a 3d vector or Tuple{2,4} for a 2 rows and 4 columns tensor,
- `T` is the type of the scalar components, should be subtype of `Number`,
- `N` is the order of the tensor, the length of `S`,
- `L` is the number of components stored internally.

`MultiValue`s are immutable.
"""
abstract type MultiValue{S,T,N,L} <: Number end

@inline Base.Tuple(arg::MultiValue) = arg.data

# Custom type printing

function show(io::IO,v::MultiValue)
  print(io,v.data)
end

function show(io::IO,::MIME"text/plain",v:: MultiValue)
  print(io,typeof(v))
  print(io,v.data)
end

###############################################################
# Introspection
###############################################################

eltype(::Type{<:MultiValue{S,T}}) where {S,T} = T
eltype(::MultiValue{S,T}) where {S,T} = T

length(::Type{<:MultiValue{S}}) where S = prod(size(MultiValue{S}))
length(a::MultiValue)  = length(typeof(a))

size(::Type{<:MultiValue{S}}) where S<:Tuple = tuple(S.parameters...)
size(a::MultiValue) = size(typeof(a))
function size(::Type{<:MultiValue{S,T,N}}, d::Integer) where {S,T,N}
    d < 1 && error("arraysize: dimension out of range")
  return d > N ? 1 :  @inbounds size(MultiValue{S})[d] # @inbounds
end
size(a::MultiValue, d::Integer) = size(typeof(a), d)

"""
    num_components(::Type{<:Number})
    num_components(a::Number)

Total number of components of a `Number` or `MultiValue`, that is 1 for scalars
and the product of the size dimensions for a `MultiValue`. This is the same as `length`.
"""
num_components(a::Number) = num_components(typeof(a))
num_components(::Type{<:Number}) = 1
num_components(T::Type{<:MultiValue}) = @unreachable "$T type is too abstract to count its components, the shape (firt parametric argument) is neede."
num_components(::Type{<:MultiValue{S}}) where S = length(MultiValue{S})

"""
    num_indep_components(::Type{<:Number})
    num_indep_components(a::Number)

Number of independent components of a `Number`, that is `num_components`
minus the number of components determined from others by symmetries or constraints.

For example, a `TensorValue{3,3}` has 9 independent components, a `SymTensorValue{3}`
has 6 and a `SymTracelessTensorValue{3}` has 5. But they all have 9 (non independent) components.
"""
num_indep_components(::Type{T}) where T<:Number = num_components(T)
num_indep_components(::T) where T<:Number = num_indep_components(T)

"""
!!! warning
    Deprecated in favor on [`num_components`](@ref).
"""
function n_components end
@deprecate n_components num_components


#######################################################################
# Other constructors and conversions implemented for more generic types
#######################################################################

zero(::Type{V}) where V<:MultiValue{S,T} where {S,T} = V(tfill(zero(T),Val(num_indep_components(V))))
zero(a::MultiValue) = zero(typeof(a))

one(::Type{V}) where V<:MultiValue = @unreachable "`one` not defined for $V"
one(a::MultiValue) = one(typeof(a))

function rand(rng::AbstractRNG,::Random.SamplerType{V}) where V<:MultiValue{D,T} where {D,T}
  Li = num_indep_components(V)
  vrand = rand(rng, SVector{Li,T})
  V(Tuple(vrand))
end


## ATM it is not possible to implement array like axes because lazy_mapping
## operations / broadcast rely on axes(::MultiValue) adopting the Number convention to return ().
#axes(::Type{<:MultiValue{S}}) where S = map(SOneTo, tuple(S.parameters...))
#axes(a::MultiValue) = axes(typeof(a))
#axes(::Type{<:MultiValue{S,T,N}},d) where {S,T,N} = d::Integer <= N ? axes(MultiValue{S})[d] : SOneTo(1)
#axes(a::MultiValue,d::Integer) = axes(typeof(a),d)

"""
    change_eltype(m::Number,::Type{T2})
    change_eltype(M::Type{<:Number},::Type{T2})

For multivalues, returns `M` or `typeof(m)` but with the component type (`MultiValue`'s parametric type `T`) changed to `T2`.

For scalars (or any non MultiValue number), `change_eltype` returns T2.
"""
change_eltype(::Type{<:Number}, T::Type) = T
change_eltype(::T,::Type{T2}) where {T<:Number,T2} = change_eltype(T,T2)
change_eltype(a::MultiValue,::Type{T2}) where T2 = change_eltype(typeof(a),T2)

get_array(a::MultiValue{S,T}) where {S,T} = convert(SArray{S,T},a)

"""
    MultiValue(a::SArray)
    MultiValue(a::MArray)

If possible (`a` needs to be of order 1, 2 or 3), converts `a` to a value of
type `<:MultiValue`.
"""
MultiValue(::StaticArray) = @notimplemented "The given StaticArray cannot be converted to a ::MultiValue."

"""
    SVector(a::MultiValue)
    SMatrix(a::MultiValue)
    SArray( a::MultiValue)
    (::Type{SA})( a::MultiValue) where SA<:StaticArray

If possible, converts `a` to a value of specified `SArray` type.

Size and element type can be specified and must match that of `a`,
e.g. `SVector{2,Int}(VectorValue(1,2))`.
"""
(::Type{SA})(a::MultiValue) where SA<:StaticArray = SA(get_array(a))

"""
    Vector(a::MultiValue)
    Matrix(a::MultiValue)
    Array( a::MultiValue)
    (::Type{A})( a::MultiValue) where A<:Array

If possible, converts `a` to a value of specified `Array` type.
The element type and number of dimension can be specified, but must be
compatible (component conversion is OK).
"""
(::Type{A})(a::MultiValue) where A<:Array = A(get_array(a))

"""
    Mutable(T::Type{<:MultiValue}) -> ::Type{<:MArray}
    Mutable(a::MultiValue)

Return the concrete mutable `MArray` type (defined by `StaticArrays.jl`) corresponding
to the `MultiValue` type T or array size and type of `a`.

See also [`mutable`](@ref).
"""
Mutable(::Type{<:MultiValue}) = @abstractmethod
# typeof(zero(...)) is for the N and L type parameters to be computed, and get
# e.g. MVector{D,T} == MMarray{Tuple{D},T,1,D} instead of MMarray{Tuple{D},T}
Mutable(::Type{<:MultiValue{S,T}}) where {S,T} = typeof(zero(MArray{S,T}))
Mutable(a::MultiValue) = Mutable(typeof(a))

"""
    mutable(a::MultiValue)

Converts `a` into a mutable array of type `MArray` defined by `StaticArrays.jl`.

See also [`Mutable`](@ref).
"""
mutable(a::MultiValue{S}) where S = MArray{S}(Tuple(get_array(a)))

###############################################################
# Conversions
###############################################################

# Direct conversion
convert(::Type{V}, arg::AbstractArray) where V<:MultiValue{S,T} where {S,T} = V(arg)
convert(::Type{V}, arg::Tuple) where V<:MultiValue{S,T} where {S,T} = V(arg)

# Inverse conversion
convert(::Type{<:NTuple{L,T}}, arg::MultiValue) where {L,T} = NTuple{L,T}(Tuple(arg))

###############################################################
# Indexing independant components
###############################################################

# This should probably not be exported, as (accessing) the data field of
# MultiValue is not a public api
"""
Previously used to transform Cartesian indices to linear indices that index `MultiValue`'s private internal storage.

!!! warning
    Deprecated, not all components of all tensors are stored anymore, so this
    index is ill defined. Use `getindex` or [`indep_comp_getindex`](@ref)
    instead of this.
"""
function data_index(::Type{<:MultiValue},i...)
  @abstractmethod
end
@deprecate data_index getindex

# The order of export of components is that of their position in the .data
# field, but the actual method "choosing" the export order is
# Gridap.Visualization._prepare_data(::Multivalue).
"""
    indep_comp_getindex(a::Number,i)

Get the `i`th independent component of `a`. It only differs from `getindex(a,i)`
when the components of `a` are interdependent, see [`num_indep_components`](@ref).
`i` should be in `1:num_indep_components(a)`.
"""
@propagate_inbounds function indep_comp_getindex(a::Number,i)
  @boundscheck @check 1 <= i <= num_indep_components(Number)
  @inbounds a[i]
end

@propagate_inbounds function indep_comp_getindex(a::T,i) where {T<:MultiValue}
  @boundscheck @check 1 <= i <= num_indep_components(T)
  @inbounds _get_data(a,i)
end

# abstraction of Multivalue data access in case subtypes of MultiValue don't
# store its data in a data field
@propagate_inbounds function _get_data(a::MultiValue,i)
  a.data[i]
end

"""
    indep_components_names(::MultiValue)

Return an array of strings containing the component labels in the order they
are exported in VTK file.

If all dimensions of the tensor shape S are smaller than 3, the components
are named with letters "X","Y" and "Z" similarly to the automatic naming
of Paraview. Else, if max(S)>3, they are labeled by integers starting from "1".
"""
function indep_components_names(::Type{MultiValue{S,T,N,L}}) where {S,T,N,L}
  return ["$i" for i in 1:L]
end

"""
    component_basis(V::Type{<:MultiValue}) -> V[ Vᵢ... ]
    component_basis(T::Type{<:Real}) -> [ one(T) ]
    component_basis(a::T<:Number)

Given a `Number` type `V` with N independent components, return a vector of
N values ``\\{ Vᵢ=V(eᵢ) \\}_i`` forming the component basis of ``\\{ u : u\\text{ isa }V\\}``,
where ``\\{eᵢ\\}_i`` is the Cartesian basis of (`eltype(V)`)ᴺ.

The `Vᵢ` verify the property that for any `u::V`,

    u = sum( indep_comp_getindex(u,i)*Vᵢ for i ∈ 1:N )
"""
component_basis(a::Number) = component_basis(typeof(a))
component_basis(T::Type{<:Number}) = [ one(T) ]
function component_basis(V::Type{<:MultiValue})
  T = eltype(V)
  Li = num_indep_components(V)
  return V[ ntuple(i -> T(i == j), Li) for j in 1:Li ]
end

"""
    representatives_of_componentbasis_dual(V::Type{<:MultiValue}) -> V[ Vᵢ... ]
    representatives_of_componentbasis_dual(T::Type{<:Real}) -> [ one(T) ]
    representatives_of_componentbasis_dual(a::V<:Number)

Given a `Number` type `V` with N independent components, return a vector of
N tensors ``\\{ Vⁱ \\}_i`` that define the form basis ``\\{ Lⁱ := (u → u ⊙ Vⁱ) \\}_i`` that
is the dual of the component basis ``\\{ Vᵢ=V(eᵢ) \\}_i`` (where ``\\{eᵢ\\}_i`` is the
Cartesian basis of (`eltype(V)`)ᴺ).

The `Vᵢ`/`Vʲ` verify the kronecker delta property ``Vᵢ⊙Vʲ = δᵢʲ``, and for any `u::V`,

    u = V( [ Lⁱ(u) for i ∈ 1:N ]... )
      = V( [ u⊙Vⁱ  for i ∈ 1:N ]... )

Rq, when `V` has dependent components, the `Vⁱ` are *not* a component basis
because ``Vʲ ≠ Vᵢ``, ``Vᵢ⊙Vⱼ≠δᵢʲ`` and

    u ≠ sum( indep_comp_getindex(u,i)*Vⁱ for i ∈ 1:N )

See also [`representatives_of_basis_dual`](@ref ) to compute a dual subspace basis.
"""
representatives_of_componentbasis_dual(a::Number) = representatives_of_componentbasis_dual(typeof(a))
representatives_of_componentbasis_dual(T::Type{<:Real}) = [ one(T) ]
function representatives_of_componentbasis_dual(V::Type{<:MultiValue})
  V = typeof(zero(V)) # makes V concrete for type inference
  N = num_indep_components(V)
  T = eltype(V)
  B = component_basis(V)

  M = MMatrix{N,N,T}(undef)
  for ci in CartesianIndices(M)
    M[ci] = B[ci[1]] ⊙ B[ci[2]]
  end
  Minv = inv(M)

  return V[ Tuple(Minv[i,:]) for i in 1:N ]
end

"""
    representatives_of_basis_dual(basis)

Same as [`representatives_of_componentbasis_dual`](@ref), but takes a basis as
arguement, i.e. a collection of basis *vectors* -- in the linear algebra sense --
of type `V<:Number`, so `basis` spans a subpace of Span``\\{ component_basis(V)...\\}``.

Computing the dual basis to a subspace basis is usefull for e.g. value types `G`
resulting from applying a derivative to function of value type `Vd` having
dependent components:

```jldoctest
Vd = SymTensorValue{2,Float64} # has 3 independent components
P = Point{2,Float64}
G = Fields.gradient_type(Vd, zero(P)) # ThirdOrderTensorValue{2, 2, 2, Float64}

# Span{ component_basis(G)... } is of dimension 8, but the gradient of a Vd
# valued function lives in the 6-dimensional space of basis:
grad_Vd_basis = [ pi⊗vi for pi in component_basis(P) for vi in component_basis(Vd)]

# of dual basis { v -> v ⊙ Gdᵢ}ᵢ for Gdᵢ in:
grad_Vd_dual_basis = representatives_of_basis_dual(grad_Vd_basis)

# kronecker delta property
[ gi ⊙ gdi for gi in grad_Vd_basis, gdi in grad_Vd_dual_basis] ≈ LinearAlgebra.I # true
```
"""
function representatives_of_basis_dual(B)
  isempty(B) && return B # also ensures we get a collection
  representatives_of_basis_dual(B...)
end
function representatives_of_basis_dual(B::T...) where {T<:Real}
  @assert length(B) = 1 && !iszero(B[1])
  T[ 1/B[1] ]
end
function representatives_of_basis_dual(B::V...) where {V<:MultiValue}
  N = length(B)
  T = eltype(V)

  M = MMatrix{N,N,T}(undef)
  for ci in CartesianIndices(M)
    M[ci] = B[ci[1]] ⊙ B[ci[2]]
  end
  Minv = inv(M)

  return V[ sum(Minv[i,:].*B) for i in 1:N ]
end

@inline function ForwardDiff.value(x::MultiValue{S,<:ForwardDiff.Dual}) where S
  VT = change_eltype(x,ForwardDiff.valtype(eltype(x)))
  data = map(ForwardDiff.value,x.data)
  return VT(data)
end
