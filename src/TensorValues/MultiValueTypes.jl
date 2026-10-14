###############################################################
# Tensor dimensions with variance
###############################################################

"""
    abstract type VariancedDimension <: Integer


"""
abstract type VariancedDimension <: Integer end

"""
    struct Cov <: VariancedDimension

`Int` wrapper representing a covariant tensor dimension. The constants `cov0 =
Cov(0)`, ..., `cov9 = Cov(9)` are exported for convenience.

`==` tests equality of the wrapped `Int` only
```julia
1 == cov1 # true
1 === cov1 # false
cov1 === Cov(1) # true
```
The wrapped `Int` can be accessed via getindex
```julia
1 === cov1[] # true
```
"""
struct Cov <: VariancedDimension
  index::Int
  Cov(index::Int) = new(index)
  Cov(index::Cov) = new(index[])
end
"""
    struct Con <: VariancedDimension

`Int` wrapper representing a contravariant tensor dimension. The constants
`con0 = Con(0)`, ..., `con9 = Con(9)` are exported for convenience.

Behaves as [`Cov`](@ref).
"""
struct Con <: VariancedDimension
  index::Int
  Con(index::Int) = new(index)
  Con(index::Con) = new(index[])
end

for i in 0:9 @eval begin
    const $(Symbol(:cov, i)) = Cov($i)
    const $(Symbol(:con, i)) = Con($i)
end end

Base.show(io::IO, ::MIME"text/plain", covi::Cov) = print(io, "cov", covi[])
Base.show(io::IO, ::MIME"text/plain", covi::Con) = print(io, "con", covi[])

# Int behavior for arithmetic of dimension (tensor size computations and checks)

Base.getindex(covi::VariancedDimension) = covi.index
Base.Int(covi::VariancedDimension) = covi[]

Base.convert(::Type{Cov}, i::Int) = Cov(i)
Base.promote_rule(::Type{Cov}, ::Type{Int}) = Int
Base.convert(::Type{Con}, i::Int) = Con(i)
Base.promote_rule(::Type{Con}, ::Type{Int}) = Int
Base.promote_rule(::Type{Con}, ::Type{Cov}) = Int

Base.:(==)(x::VariancedDimension, y::Integer) = x[] == y
Base.:(==)(x::Integer, y::VariancedDimension) = x == y[]
Base.:(==)(x::VariancedDimension, y::VariancedDimension) = x[] == y[]
Base.:(<)(x::VariancedDimension, y::Integer) = x[] < y
Base.:(<)(x::Integer, y::VariancedDimension) = x < y[]
Base.:(<)(x::VariancedDimension, y::VariancedDimension) = x[] < y[]
Base.:(>>)(x::VariancedDimension, y::Int) = x[] >> y
Base.trailing_zeros(x::VariancedDimension) = trailing_zeros(x[])

Base.:+(a::Cov, b::Cov) = Cov(a[] + b[])
Base.:-(a::Cov, b::Cov) = Cov(a[] - b[])
Base.:*(a::Cov, b::Cov) = Cov(a[] * b[])
Base.:/(a::Cov, b::Cov) = Cov(a[] / b[])
Base.div(a::Cov, b::Cov, r::RoundingMode=RoundToZero) = div(a[], b[], r)
Base.:-(x::Cov) = Cov(-x[])

Base.:+(a::Con, b::Con) = Con(a[] + b[])
Base.:-(a::Con, b::Con) = Con(a[] - b[])
Base.:*(a::Con, b::Con) = Con(a[] * b[])
Base.:/(a::Con, b::Con) = Con(a[] / b[])
Base.div(a::Con, b::Con, r::RoundingMode=RoundToZero) = div(a[], b[], r)
Base.:-(x::Con) = Con(-x[])

"""
    Dim = Union{Int, Cov, Con}

Valid type of values to use in `MultiValue`'s size, `S=Tuple{I...}` where each `D` in `I` is a `Dim`.
"""
const Dim = Union{Int, Cov, Con}

# convert a size Tuple S to the tuple size(::MultiValue{S})
# @generated because the functions needt to read the `parameters` field of the
# S::DataType, while remaining constant folded
@generated function _Size_to_size(::Type{S}) where S<:Tuple
  s = tuple(S.parameters...)
  isempty(filter(D -> !isa(D, Dim), s)) || return :(
    throw(ErrorException("$S is not a tensor size Tuple, it should only contain [`Dim`](@ref) indices."))
  )
  :( $s )
end

"""
    const _Scalar = Union{Real,Complex}

Abstract type for the scalar types that `MultiValue` support component-wise
operations with.
"""
const _Scalar = Union{Real,Complex}

"""
    is_same_variance(i::Dim, j::Dim) = typeof(i) === typeof(j)


Wether dimensions `i` and `j` have the same variance.
"""
is_same_variance(i::Dim, j::Dim) = typeof(i) === typeof(j)

"""
    is_opposit_variance(i::Dim, j::Dim)


Wether dimensions `i` and `j` have variance compatible for contraction, both
`Int` or `Cov`↔`Con`.

This is not the contrary of `is_same_variance`, because e.g.
    is_opposite_variance(::Int, ::Con) == false
"""
is_opposit_variance(i::Int, j::Int) = true
is_opposit_variance(i::Con, j::Cov) = true
is_opposit_variance(i::Cov, j::Con) = true
is_opposit_variance(i::Dim, j::Dim) = false

"""
    can_contract(i::Dim, j::Dim) = is_opposit_variance(i,j) && i==j

Whether dimensions `i` and `j` can be contracted.
"""
can_contract(i::Dim, j::Dim) = is_opposit_variance(i,j) && i==j

@inline function _ensure_contract(i::Dim, j::Dim)
  can_contract(i, j) || throw(DimensionMismatch(
    "Indices $i and $j cannot be contracted, see `can_contract`."))
  nothing
end

"""
    inverse_variance(i::Dim) -> Dim
    inverse_variance(s::Tuple) -> tuple
    inverse_variance(::Type{S<:Tuple}) -> Tuple

Dimension `i` with its variance flipped, `Cov`↔`Con` or `Int`↔`Int`. If a size
tuple or Size Tuple is given, `inverse_variance` is applied to each index.

`can_contract(i, inverse_variance(i))` holds for every `i::Dim`.
"""
inverse_variance(i::Int) = i
inverse_variance(i::Cov) = Con(i[])
inverse_variance(i::Con) = Cov(i[])
inverse_variance(s::Tuple) = map(inverse_variance, s)
inverse_variance(::Type{S}) where S<:Tuple = Tuple{map(inverse_variance, _Size_to_size(S))...}

_rm_size_variance(s::Tuple) = tuple(Int.(s)...)
_rm_size_variance(::Type{S}) where S<:Tuple = Tuple{Int.(_Size_to_size(S))...}


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
eltype(V::Type{<:MultiValue}) = @unreachable "$V type has no element type T."

length(::Type{<:MultiValue{S}}) where S = prod(_rm_size_variance(size(MultiValue{S})))
length(a::MultiValue) = length(typeof(a))
length(V::Type{<:MultiValue}) = @unreachable "$V type is too abstract to count its components, the size information is needed."

size(::Type{<:MultiValue{S}}) where S<:Tuple = _Size_to_size(S)
size(a::MultiValue) = size(typeof(a))
function size(::Type{<:MultiValue{S,T,N}}, d::Integer) where {S,T,N}
    d < 1 && error("arraysize: dimension out of range")
  return d > N ? 1 :  @inbounds size(MultiValue{S})[d] # @inbounds
end
size(a::MultiValue, d::Integer) = size(typeof(a), d)
size(V::Type{<:MultiValue}) = @unreachable "$V type has no size information."

"""
    tensor_type(a::Number)
    tensor_type(::Type{V<:Number})

The triple `(r, s, u)` where
- `s` the number of covariant dimensions, i.e. indices that `isa Cov`,
- `r` the number of contravariant dimensions, i.e. indices that `isa Con`,
- `u` the number of un-varianced dimensions, i.e. indices that `isa Int`.
"""
tensor_type(a::V)      where V<:Number     = _tensor_type(size(a))
tensor_type(::Type{V}) where V<:MultiValue = _tensor_type(size(V))
function _tensor_type(sizetuple::S) where S<:Tuple
  @assert all(i -> isa(i,Dim), sizetuple) "Invalid dimension indices in size: $sizetuple"
  r = count(i -> isa(i,Con), sizetuple)
  s = count(i -> isa(i,Cov), sizetuple)
  u = length(sizetuple) - s - r
  (r,s,u)
end

"""
    num_components(::Type{<:Number})
    num_components(a::Number)

Total number of components of a `Number` or `MultiValue`, that is 1 for scalars
and the product of the size dimensions for a `MultiValue`. This is the same as `length`.
"""
num_components(a::Number) = num_components(typeof(a))
num_components(::Type{<:Number}) = 1
num_components(V::Type{<:MultiValue}) = @unreachable "$V type is too abstract to count its components, the shape (first parametric argument) is needed."
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
    get_indep_components(a::Number)

Tuple of the [`num_indep_components`](@ref) independent components of `a`, in the
same order as [`indep_comp_getindex`](@ref).
"""
get_indep_components(a::Number) = (a,)
get_indep_components(a::MultiValue) = Tuple(a)

"""
    num_stored_components(::Type{<:MultiValue})
    num_stored_components(a::MultiValue{S,T,N,L}) = L

Return the number of stored components of a tensor.

An example for which this differs from `num_indep_components` is
`SymTracelessTensorValue`, that stores the last diagonal value although it
equals minus the sum of the others.
"""
num_stored_components(V::Type{<:MultiValue}) = num_indep_components(V)
num_stored_components(a::MultiValue) = num_stored_components(typeof(a))

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
zero(V::Type{<:MultiValue}) = @unreachable "$V type is too abstract, the size and eltype parameters."

one(::Type{V}) where V<:MultiValue = @unreachable "`one` not defined for $V"
one(a::MultiValue) = one(typeof(a))

function rand(rng::AbstractRNG,::Random.SamplerType{V}) where V<:MultiValue{D,T} where {D,T}
  Li = num_indep_components(V)
  vrand = rand(rng, SVector{Li,T})
  V(Tuple(vrand))
end

"""
    make_concretetype(::Type{T}) where T <: Number

Return a concrete number type for `T`, that is `T` itself if it is already
concrete, and `typeof(zero(T))` otherwise.
"""
function make_concretetype(::Type{T}) where T <: Number
  TT = ifelse(isconcretetype(T),T,typeof(zero(T)))
  @check isconcretetype(TT) "Type $(T) cannot be made concrete."
  return TT
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

"""
    change_variance(m::MultiValue, ::Val{Sv})
    change_variance(M::Type{<:MultiValue}, ::Val{Sv})
    change_variance(m::MultiValue, ::Type{S})
    change_variance(M::Type{<:MultiValue}, ::Type{S})

Return `M` or `typeof(m)` with its dimensions replaced by those of the tuple
`Sv`, or of the size `Tuple` `S`. The returned type is not concrete, but can
un-ambiguously be made so using [`make_concretetype`](@ref).
"""
change_variance(a::MultiValue, S) = change_variance(typeof(a), S)
change_variance(V::Type{<:MultiValue}, ::Type{S}) where S<:Tuple =
  change_variance(V, Val(_Size_to_size(S)))

function _check_change_variance(::Type{V}, Sv::Tuple, uniform::Bool=false) where V<:MultiValue
  sz = size(V)
  length(Sv) == length(sz) || throw(ArgumentError("$V has $(length(sz)) indices,
    its variance cannot be changed to the $(length(Sv)) dimensions $Sv."))
  all(k -> Int(Sv[k]) == Int(sz[k]), eachindex(Sv)) || throw(ArgumentError(
    "`change_variance` only flips variances, but $Sv and the size $sz of $V do
     not have the same dimensions."))
  uniform && !isempty(Sv) && !all(d -> d === first(Sv), Sv) && throw(ArgumentError(
    "$V carries a single dimension for several of its indices, so no type
     expresses it with the dimensions $Sv."))
  nothing
end

get_array(a::MultiValue{S,T}) where {S,T} = convert(SArray{_rm_size_variance(S),T},a)

"""
    MultiValue(a::SArray)
    MultiValue(a::MArray)

Converts `a` to a tensor of type `<:MultiValue`.
""" # the S param is there to allow NTensorValue's specialization to be more specific
MultiValue(::StaticArray{S}) where S = @abstractmethod

"""
    SVector(a::MultiValue)
    SMatrix(a::MultiValue)
    SArray( a::MultiValue)
    (::Type{SA})( a::MultiValue) where SA<:StaticArray

If possible, converts `a` to a value of specified `SArray` type.

Size and element type can be specified and must match that of `a`,
e.g. `SVector{2,Int}(VectorValue(1,2))`.

For `MultiValue` types with dependent component, the conversion works but the
dependency information is lost.
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
Mutable(::Type{<:MultiValue{S,T}}) where {S,T} = typeof(zero(MArray{_rm_size_variance(S),T}))
Mutable(a::MultiValue) = Mutable(typeof(a))

"""
    mutable(a::MultiValue)

Converts `a` into a mutable array of type `MArray` defined by `StaticArrays.jl`.

See also [`Mutable`](@ref).
"""
mutable(a::MultiValue{S}) where S = MArray{_rm_size_variance(S)}(Tuple(get_array(a)))

###############################################################
# Conversions
###############################################################

# Direct conversion
convert(::Type{V}, arg::AbstractArray) where V<:MultiValue{S,T} where {S,T} = V(arg)
convert(::Type{V}, arg::Tuple) where V<:MultiValue{S,T} where {S,T} = V(arg)

# Inverse conversion
convert(::Type{<:NTuple{L,T}}, arg::MultiValue) where {L,T} = NTuple{L,T}(Tuple(arg))

function reinterpret(a::Array{V}) where V<:MultiValue{S,T,N,L} where {S,T,N,L}
  b = reinterpret(T,a)
  sa = size(a)
  s = (L,sa...)
  reshape(b,s)
end

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

# similar to CartesianIndices, but indicating the linear index in the stored
# component, with the sign flip to apply (for e.g. skew tensors). `nothing` for a
# type that materialises a component no single full tensor index reads, which has
# no such map.
function stored_comp_indices(V::Type{<:MultiValue})
  Li = num_indep_components(V)
  Li == num_stored_components(V) || @notimplemented "$V stores components that
    are not independent, it needs a `stored_comp_indices` method for this operation."
  V(Tuple(1:Li))
end

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
  @boundscheck @check 1 <= i <= num_indep_components(typeof(a))
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
    indep_components_names(::Type{<:MultiValue{S}})

Return an array of strings containing the component labels in the order they
are exported in VTK file.

If all dimensions of the tensor shape `S` are smaller than 3, the components
are named with letters `"X"`,`"Y"` and `"Z"` similarly to the automatic naming
of Paraview. Else, if `maximum(S)>3`, they are labeled by integers starting from `"1"`.
"""
indep_components_names(a::MultiValue) = indep_components_names(typeof(a))
function indep_components_names(::Type{V}) where V<:MultiValue{S,T,N,L} where {S,T,N,L}
  s = size(NTensorValue{S})

  # if there are dependant components, numbered linearly by default
  if L != prod(s)
    @warn "indep_components_names has not been implemented for a tensor type
    with dependent components, falling back to default names [string(i) for i in 1:L]"
    return [string(i) for i in 1:L]
  end

  indep_components_names(MultiValue{S})
end

function indep_components_names(::Type{V}) where V<:MultiValue{S} where {S}
  s = size(NTensorValue{S})

  if maximum(s; init=0) > 3
    names = [ join(Tuple(ci)) for ci in CartesianIndices(s) ]
  else
    c_name = ["X", "Y", "Z"]
    c_names(ci) = join(c_name[i] for i in Tuple(ci))
    names =  [ c_names(ci) for ci in CartesianIndices(s) ]
  end
  reshape(names, length(names))
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
    representatives_of_componentbasis_dual(V::Type{<:MultiValue}, metric_args...) -> V[ Vᵢ... ]
    representatives_of_componentbasis_dual(T::Type{<:Real}, metric_args...) -> [ one(T) ]
    representatives_of_componentbasis_dual(a::V<:Number, metric_args...)

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

`metric_args` is the metric information the inner product of a varianced `V`
requires, see [`parse_metric_data`](@ref) and [`inner`](@ref).

See also [`representatives_of_basis_dual`](@ref ) to compute a dual subspace basis.
"""
representatives_of_componentbasis_dual(a::Number, metric_args...) =
  representatives_of_componentbasis_dual(typeof(a), metric_args...)
representatives_of_componentbasis_dual(T::Type{<:Real}, metric_args...) = [ one(T) ]
function representatives_of_componentbasis_dual(V::Type{<:MultiValue}, metric_args...)
  V = typeof(zero(V)) # makes V concrete for type inference
  N = num_indep_components(V)
  iszero(N) && return V[]
  Minv = inv(_gram_matrix(component_basis(V), metric_args...))

  return V[ Tuple(Minv[i,:]) for i in 1:N ]
end

# The Gram matrix `Mᵢⱼ = Bᵢ⊙Bⱼ` of a non-empty family of tensors, the inner
# product being taken with the given metric.
function _gram_matrix(B, metric_args...)
  N = length(B)
  m11 = inner(first(B), first(B), metric_args...)
  M = MMatrix{N,N,typeof(m11)}(undef)
  for ci in CartesianIndices(M)
    M[ci] = inner(B[ci[1]], B[ci[2]], metric_args...)
  end
  M
end

"""
    representatives_of_basis_dual(basis, metric_args...)

Same as [`representatives_of_componentbasis_dual`](@ref), but takes a basis as
arguement, i.e. a collection of basis *vectors* -- in the linear algebra sense --
of type `V<:Number`, so `basis` spans a subpace of Span{`component_basis(V)...`}.

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
function representatives_of_basis_dual(B, metric_args...)
  isempty(B) && return B
  _representatives_of_basis_dual(B, metric_args...)
end

function _representatives_of_basis_dual(B, metric_args...)
  V = eltype(B)
  N = length(B)
  Minv = inv(_gram_matrix(B, metric_args...))

  return V[ sum(Minv[i,:].*B) for i in 1:N ]
end
function _representatives_of_basis_dual(B::AbstractVector{T}) where {T<:Real}
  @assert isone(length(B)) && !iszero(first(B))
  T[ 1/first(B) ]
end

@inline function ForwardDiff.value(x::MultiValue{S,<:ForwardDiff.Dual}) where S
  VT = change_eltype(x,ForwardDiff.valtype(eltype(x)))
  data = map(ForwardDiff.value,x.data)
  return VT(data)
end

###############################################################
# Combinations
###############################################################

"""
    sorted_combinations(D,k; right_to_left=false)

Return a (collected) vector of all the combinations I_i of {1:`D`} of length `k`:

    1 ≤ I\\_1 < ... < I\\_k ≤ `D`

sorted in (left-to-right, i.e. standard) lexicographic order, e.g.

```julia
[ [1,2], [1,3], [2,3] ]  # for D=3, k=2\\
[ [1,2], [1,3], [1,4], [2,3], [2,4], [3,4] ]  # for D=4, k=2
```

If `right_to_left` is `true`, the combinations are compared from their last
index to their first instead, e.g.

```julia
[ [1,2], [1,3], [2,3], [1,4], [2,4], [3,4] ]  # for D=4, k=2
```

This is the order in which `get_faces` numbers the sub-faces of a simplex.

`k=0` gives the single empty combination `[ Int[] ]`, while `k<0` or `k>D` give
no combination at all.

See also [`combination_index`](@ref).
"""
function sorted_combinations(D::Int,k::Int; right_to_left=false)
  combis = collect(combinations(1:D,k)) # empty if k < 0 or k > D
  right_to_left && sort!(combis, by=reverse)
  return combis
end

# Multiplicative form of `binomial`: `c * (n-k+i) ÷ i`, unlike `binomial` it
# inlines and carries no overflow check (which are irrelevant for tensor type).
@inline function _unchecked_binomial(n, k)
  c = 1
  for i in 1:k; c = c * (n-k+i) ÷ i; end
  c
end

"""
    combination_index(I, D; right_to_left=false)

Linear index of `I` amongst the combinations of 1:`D` of the same size `k`,
sorted like in [`sorted_combinations`](@ref), that is

    (combination_index(I, D; right_to_left), I) ∈ enumerate(sorted_combinations(D,k; right_to_left))

The right-to-left index does not depend on `D`, unlike the default left-to-right one.
"""
@inline function combination_index(combi, D; right_to_left=false)
  @check issorted(combi)
  @check isempty(combi) || last(combi) ≤ D
  k = length(combi)
  # Counting the combinations that come before `combi`: those first differing
  # from it at position i, read from the right, number binomial(combi[i]-1, i).
  right_to_left && return sum(_unchecked_binomial(combi[i]-1, i) for i in 1:k; init=0) + 1
  # Counting the combinations that come after `combi`: those first exceeding it
  # at position i number binomial(D-combi[i], k-i+1).
  return _unchecked_binomial(D,k) - sum(_unchecked_binomial(D-combi[i], k-i+1) for i in 1:k; init=0)
end

"""
    combination_complement(I, D)

Given a `k`-combination `I` of elements of `1:D`, return the unique (`D`-`k`)
combination of `1:D` "`Icomp`" such that `I ∪ Icomp ⊇ 1:D`.

See also [`combination_sign`](@ref).
"""
function combination_complement(combi, D)
  k = length(combi)
  combi_comp = Vector{Int}(undef, D-k)
  curr_combi, curr_comp = 1, 1
  for i in 1:D
    if curr_combi ≤ k && combi[curr_combi] == i
      curr_combi += 1
    else
      combi_comp[curr_comp] = i
      curr_comp += 1
    end
  end
  combi_comp
end

"""
    combination_sign(I)

Given a combination `I`, return the sign of the permutation resulting from the
concatenation of `I` and its complement [`combination_complement(I,D)`](@ref
combination_complement).
"""
function combination_sign(combi)
  @check issorted(combi)
  i, k, acc, delta = 1, 1, 0, 0
  while k <= length(combi)
    if combi[k] == i
      acc += delta
      k += 1
    else
      delta += 1
    end
    i += 1
  end
  return iseven(acc) ? 1 : -1
end

"""
    sub_combinations_ids(J, D)

Return a vector containing the [`combination_index`](@ref) of the `k` different
combinations `J\\J[i]` of 1:`D` for 1 ≤ i ≤ `k`, where `k=length(J)`.

`D` is the ambient dimension the indices are taken in, that is the one `J`
itself is a combination of — not the dimension of the space the forms live in.

See also [`sup_combinations_ids`](@ref), the inverse direction.
"""
function sub_combinations_ids(combi, D)
  k = length(combi)
  sub_combi = Vector{Int}(undef, k-1)
  sub_combi_ids = Vector{Int}(undef, k)
  for i in 1:k
    for j in 1:k-1
      sub_combi[j] = combi[j + Int(j≥i)]
    end
    sub_combi_ids[i] = combination_index(sub_combi, D)
  end
  sub_combi_ids
end

"""
    sup_combinations_ids(I, D)

Return a vector containing the [`combination_index`](@ref) of the `D`-`k`
different combinations `I ∪ {j}` of 1:`D` for `j ∉ I`, where `k=length(I)`,
ordered by increasing `j`, which is also the lexicographic order of the
combinations.
"""
function sup_combinations_ids(combi, D)
  k = length(combi)
  sup_combi = Vector{Int}(undef, k+1)
  sup_combi[2:end] .= combi
  sup_combi_ids = Vector{Int}(undef, D-k)
  p, n = 1, 1
  for j in 1:D
    sup_combi[p] = j
    if p ≤ k && combi[p] == j
      p += 1
    else
      sup_combi_ids[n] = combination_index(sup_combi, D)
      n += 1
    end
  end
  sup_combi_ids
end

"""
    sorted_multisets(D,N)

Return a (collected) vector of all the non-decreasing multi-indices I_i of {1:`D`} of length `N`:

    1 ≤ I\\_1 ≤ ... ≤ I\\_N ≤ `D`

sorted in (left-to-right, i.e. standard) lexicographic order, e.g.

```julia
[ [1,1], [1,2], [1,3], [2,2], [2,3], [3,3] ]  # for D=3, N=2
```

These index the independent components of a fully symmetric `N`th-order tensor,
of which there are `binomial(D+N-1,N)`.

`N=0` gives the single empty multi-index `[ Int[] ]`, while `N<0`, or `D=0` with
`N>0`, give no multi-index at all.

See also [`multiset_index`](@ref), [`sorted_combinations`](@ref).
"""
function sorted_multisets(D::Int,N::Int)
  # A non-decreasing I of {1:D} is the strictly increasing (I_i + i - 1) of
  # {1:D+N-1}, which is the bijection `multiset_index` inverts.
  [ combi .- (0:N-1) for combi in sorted_combinations(D+N-1, N) ]
end

"""
    multiset_index(I, D)

Linear index of the non-decreasing multi-index `I` amongst the multi-indices of
1:`D` of the same length `N`, sorted like in [`sorted_multisets`](@ref), that is

    (multiset_index(I, D), I) ∈ enumerate(sorted_multisets(D,N))

This is the index of the component `a[I...]` in the independent components of a
fully symmetric tensor `a` of order `N` and dimension `D`.

See also [`combination_index`](@ref).
"""
@inline function multiset_index(multiset, D)
  @check issorted(multiset)
  @check isempty(multiset) || last(multiset) ≤ D
  N = length(multiset)
  # `combination_index` of the associated combination (I_i + i - 1) of 1:D+N-1,
  # with the shift folded into the binomials.
  return binomial(D+N-1,N) - sum(_unchecked_binomial(D+N-i-Ii, N-i+1) for (i,Ii) in enumerate(multiset); init=0)
end

"""
    sorting_sign(inds::Integer...)

Signature of the permutation sorting `inds`, or `0` if `inds` repeats an index.

This is `levicivita(sortperm(inds))` extended by `0` to non-injective `inds`.
"""
@inline function sorting_sign(inds::Vararg{Integer,N}) where N
  ninv = 0
  for j in 1:N, k in j+1:N
    inds[k] == inds[j] && return 0
    inds[k] <  inds[j] && (ninv += 1)
  end
  iseven(ninv) ? 1 : -1
end
