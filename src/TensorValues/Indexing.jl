# It would be better to keep IndexLinear style for VectorValue, TensorValue and
# ThirdOrderTensorValue, their eachindex isn't as efficient as possible
Base.IndexStyle(::MultiValue) = IndexCartesian()
Base.IndexStyle(::Type{<:MultiValue}) = IndexCartesian()

# Necessary overloads due to wrong ::Number defaults
lastindex(arg::MultiValue) = length(arg)
lastindex(arg::MultiValue, d::Integer) = (@inline; size(arg, d)) # generic
lastindex(arg::MultiValue, d::Int)     = (@inline; size(arg, d)) # method ambiguity with base

# Gridap broadcast of some operation on ::MultiValue rely on Base.axes adopting
# the Number convension (all MultiValue have axes `()` )
# This is the AbstractArray like axes
_axes(arg::MultiValue) = map(SOneTo, size(arg))

CartesianIndices(arg::MultiValue) = CartesianIndices(_axes(arg)) # Array axes
LinearIndices(arg::MultiValue) = LinearIndices(_axes(arg))

eachindex(arg::MultiValue) = CartesianIndices(arg)
eachindex(::IndexCartesian, arg::MultiValue) = eachindex(arg)
eachindex(::IndexLinear, arg::MultiValue) = SOneTo(length(arg))

keys(arg::MultiValue) = CartesianIndices(map(SOneTo, size(arg)))
keys(s::IndexStyle, arg::MultiValue) = eachindex(s, arg)

"""
    getindex(arg::MultiValue, inds...)
    getindex(arg::MultiValue, i::Integer)
    getindex(arg::MultiValue) = arg

`MultiValue`s support a large subset of the indexing methods of
`AbstractArray`s and `StaticArray`s, including `inds` argument that is a mixed
tuples of `Integer`, `CartesianIndex`, slices and common range and array types.

The `Number` convention is used when no indices are provided: `arg[]` returns `arg`.

# Examples

```julia
julia> t = TensorValue(1:9...)
TensorValue{3, 3, Int64, 9}(1, 2, 3, 4, 5, 6, 7, 8, 9)

julia> t[1,end]
7

julia> t[:,2]
VectorValue{3, Int64}(4, 5, 6)

julia> t[ isodd.(SMatrix(t))]
5-element Vector{Int64}:
 1
 3
 5
 7
 9
```

# Extended help

Similarly to StaticArray.jl, the type of the returned value depend on
inferability of the size from the type of the indices in `inds`.

Indeed, `getindex` returns
- a scalar (component) if `inds` contains only "scalars", i.e. `Integer` and `CartesianIndex`,
- a `<:MultiValue` tensor if `inds` additionally contains statically inferable index ranges/sets such as `Colon()`/`:`, `SOneTo` or `StaticArray`,
- an `Array` if `inds` additionally contains dynamic index ranges/sets such as `UnitRange`, `OneTo` or other types of `AbstractArray`.

# warning
  Indexing methods that do not return a scalars will loose any symmetry or other
  known internal constraint, as in `SymTensorValue(1:6...)[SOneTo(2),SOneTo(2)]`

```julia
julia> t[SOneTo(2),SOneTo(3)]
TensorValue{2, 3, Int64, 6}(1, 2, 4, 5, 7, 8)

julia> t[Base.OneTo(2),Base.OneTo(3)]
2×3 Matrix{Int64}:
 1  4  7
 2  5  8

julia> t[MVector(1,2),CartesianIndex(3)]
VectorValue{2, Int64}(7, 8)

julia> t[1:2,CartesianIndex(3)]
2-element Vector{Int64}:
 7
 8


julia> mask = [true false false; false true false; false false true]
3×3 Matrix{Bool}:
 1  0  0
 0  1  0
 0  0  1

julia> t[mask]
3-element Vector{Int64}:
 1
 5
 9
```

""" # Those methods are necessary because MultiValue subtypes Number instead of AbstractArray
@propagate_inbounds getindex(arg::MultiValue, i::Integer) = getindex(arg, CartesianIndices(arg)[i])
@propagate_inbounds getindex(arg::MultiValue) = arg
# Size-inferable "scalar" indexing
const _ScalarIndices = Union{Integer, CartesianIndex}
@propagate_inbounds getindex(arg::MultiValue, inds::_ScalarIndices...) = getindex(arg, to_indices(arg, inds)...)
# Method to avoid infinite recursion in case of wrong number of scalar indices
@propagate_inbounds getindex(arg::MultiValue, inds::Integer...) = (checkbounds(arg,inds...); @unreachable)
# Size-inferable "array" indexing,
const _StaticIndices = Union{Colon,SOneTo,StaticArray{<:Tuple, <:Union{Int32,Int64}},_ScalarIndices}
# the conversion to SArray is only necessary when there are CartesianIndex in `inds`, see https://github.com/JuliaArrays/StaticArrays.jl/issues/1059
@propagate_inbounds getindex(arg::MultiValue, inds::_StaticIndices...) = MultiValue(SArray(getindex(get_array(arg), inds...)))
# Not size-inferable "array" indexing, returns ::Base.Array. Includes ::OneTo, Array{Bool}, etc.
const _DynamicIndices= Union{AbstractArray,_StaticIndices}
@propagate_inbounds getindex(arg::MultiValue, inds::_DynamicIndices...) = getindex(get_array(arg), inds...)


# Cartesian indexing style implementation
@propagate_inbounds function getindex(arg::VectorValue,i::Integer)
  @boundscheck @check checkbounds(arg,i) === nothing
  @inbounds arg.data[i]
end

@propagate_inbounds function getindex(arg::TensorValue{D},i::Integer,j::Integer) where D
  @boundscheck @check checkbounds(arg,i,j) === nothing
  index = _2d_tensor_linear_index(D,i,j)
  @inbounds arg.data[index]
end

@propagate_inbounds function getindex(arg::AbstractSymTensorValue{D},i::Integer,j::Integer) where D
  @boundscheck @check checkbounds(arg,i,j) === nothing
  index = _2d_sym_tensor_linear_index(D,i,j)
  @inbounds arg.data[index]
end

@propagate_inbounds function getindex(arg::SkewSymTensorValue{D,T},i::Integer,j::Integer) where {D,T}
  @boundscheck @check checkbounds(arg,i,j) === nothing
  i == j && return zero(T)
  index = _2d_skew_sym_tensor_linear_index(D,i,j)
  v = @inbounds arg.data[index]
  i<j ? v : -v
end

@propagate_inbounds function getindex(arg::SymFourthOrderTensorValue{D},i::Integer,j::Integer,k::Integer,l::Integer) where D
  @boundscheck @check checkbounds(arg,i,j,k,l) === nothing
  index = _4d_sym_tensor_linear_index(D,i,j,k,l)
  @inbounds arg.data[index]
end

@propagate_inbounds function getindex(arg::HighOrderTensorValue{S,T,N}, inds::Vararg{Integer,N}) where {S,T,N}
  @boundscheck @check checkbounds(arg, inds...) === nothing
  @inbounds arg.data[inds...]
end
# The following attempt to switch to IndexLinear style failed due to the manual
# indexing dispatches above...
#
#Base.IndexStyle(::HighOrderTensorValue) = IndexLinear()
#Base.IndexStyle(::Type{<:HighOrderTensorValue}) = IndexLinear()
#Base.checkbounds(arg::HighOrderTensorValue, i::Integer) = checkbounds(arg.data, i)
#@propagate_inbounds function getindex(arg::HighOrderTensorValue, i::Integer)
#  @boundscheck @check checkbounds(arg, i) === nothing
#  @inbounds @inline getindex(arg.data, i)
#end


function Base.checkbounds(A::MultiValue{S}, I::Integer...) where S
  if CartesianIndex(I...) ∉ CartesianIndices(A)
    throw(BoundsError(A,I))
  end
  nothing
end

@inline iterate(arg::MultiValue)        = iterate(arg.data)
@inline iterate(arg::MultiValue, state) = iterate(arg.data, state)

# to remove next major
data_index(::Type{<:VectorValue},i) = i
data_index(::Type{<:TensorValue{D}},i,j) where D = _2d_tensor_linear_index(D,i,j)
data_index(::Type{<:AbstractSymTensorValue{D}},i,j) where D = _2d_sym_tensor_linear_index(D,i,j)
data_index(::Type{<:ThirdOrderTensorValue{D1,D2}},i,j,k) where {D1,D2} = _3d_tensor_linear_index(D1,D2,i,j,k)
data_index(::Type{<:SymFourthOrderTensorValue{D}},i,j,k,l) where D = _4d_sym_tensor_linear_index(D,i,j,k,l)

_symmetric_index_gaps(i::Integer) = i*(i-1)÷2
_skew_symetric_index_gaps(i::Integer) = i*(i+1)÷2

_2d_tensor_linear_index(D,i,j) = ((j-1)*D)+i

_3d_tensor_linear_index(D1,D2,i,j,k) = (k-1)*D1*D2+(j-1)*D1+i

function _2d_sym_tensor_linear_index(D,i,j)
  _j,_i = minmax(i,j)
  index=_2d_tensor_linear_index(D,_i,_j)-_symmetric_index_gaps(_j)
  index
end

function _2d_skew_sym_tensor_linear_index(D,i,j)
  _j,_i = minmax(i,j)
  index=_2d_tensor_linear_index(D,_i,_j)-_skew_symetric_index_gaps(_j)
  index
end

#function _4d_sym_tensor_linear_index(D,i,j,k,l)
#  _j=min(i,j)
#  _i=max(i,j)
#  _l=min(k,l)
#  _k=max(k,l)
#  block_length=_symmetric_index_gaps(D+1)
#  element_index=_2d_tensor_linear_index(D,_i,_j)-_symmetric_index_gaps(_j)
#  block_index=_2d_tensor_linear_index(D,_l,_k)-_symmetric_index_gaps(_l)
#  index=(block_index-1)*block_length+element_index
#  index
#end

function _4d_sym_tensor_linear_index(D,i,j,k,l)
  block_length = (D*(D+1))÷2
  block_index = _2d_sym_tensor_linear_index(D,i,j)
  element_index = _2d_sym_tensor_linear_index(D,k,l)
  index=(block_index-1)*block_length+element_index
  index
end

###############################################################
# Voigt and Mandel notation
###############################################################

# Compile-time helpers (called only inside @generated function bodies)
function _voigt_pairs(D)
  pairs = NTuple{2,Int}[]
  for i in 1:D; push!(pairs, (i,i)); end
  for i in 1:D, j in i+1:D; push!(pairs, (i,j)); end
  pairs
end

function _voigt_inv(D)
  pairs = _voigt_pairs(D)
  inv_map = zeros(Int, D, D)
  for (k, (i,j)) in enumerate(pairs)
    inv_map[i,j] = k; inv_map[j,i] = k
  end
  inv_map
end

# Promoted eltype for Mandel (off-diagonal entries multiplied by √2)
_mandel_eltype(::Type{T}) where {T} =
  Base.promote_op(*, typeof(sqrt(2*one(T))), T)

"""
    to_voigt(a::SymTensorValue{D,T,L}) -> VectorValue{L,T}

Encode a symmetric second-order tensor in Voigt notation. Components are
ordered with diagonals first — (1,1),(2,2),…,(D,D) — followed by
upper-triangular off-diagonal pairs (i,j) with i<j in lexicographic order.
"""
@generated function to_voigt(a::SymTensorValue{D,T,L}) where {D,T,L}
  pairs = _voigt_pairs(D)
  str = join(["a[$(p[1]),$(p[2])], " for p in pairs])
  Meta.parse("VectorValue{$L,T}(($str))")
end

"""
    from_voigt(v::VectorValue{L,T}) -> SymTensorValue

Decode a Voigt-encoded vector to a symmetric second-order tensor.
`L` must be a triangular number, i.e. L = D(D+1)/2 for some integer D.
"""
@generated function from_voigt(v::VectorValue{L,T}) where {L,T}
  D = (isqrt(1+8*L)-1) ÷ 2
  @assert L == D*(D+1)÷2 "from_voigt: VectorValue length $L is not a valid Voigt vector length"
  inv_map = _voigt_inv(D)
  str = join(["v[$(inv_map[i,j])], " for i in 1:D for j in i:D])
  Meta.parse("SymTensorValue{$D,T}(($str))")
end

"""
    to_mandel(a::SymTensorValue{D,T,L}) -> VectorValue{L,S}

Encode a symmetric second-order tensor in Mandel notation. The ordering
follows Voigt convention, but off-diagonal components are scaled by √2 so
that `to_mandel(a) ⋅ to_mandel(b) == inner(a, b)` for any `a`, `b`.
The output eltype `S` is `promote_op(*, typeof(√2), T)`.
"""
@generated function to_mandel(a::SymTensorValue{D,T,L}) where {D,T,L}
  S = _mandel_eltype(T)
  pairs = _voigt_pairs(D)
  terms = [i==j ? "convert($S, a[$i,$j])" : "sqrt($S(2)) * a[$i,$j]" for (i,j) in pairs]
  str = join(terms, ", ")
  Meta.parse("VectorValue{$L,$S}(($str,))")
end

"""
    from_mandel(v::VectorValue{L,T}) -> SymTensorValue

Decode a Mandel-encoded vector to a symmetric second-order tensor.
`L` must be a triangular number, i.e. L = D(D+1)/2 for some integer D.
"""
@generated function from_mandel(v::VectorValue{L,T}) where {L,T}
  D = (isqrt(1+8*L)-1) ÷ 2
  @assert L == D*(D+1)÷2 "from_mandel: VectorValue length $L is not a valid Mandel vector length"
  S = _mandel_eltype(T)
  inv_map = _voigt_inv(D)
  terms = [i==j ? "convert($S, v[$(inv_map[i,j])])" : "v[$(inv_map[i,j])] / sqrt($S(2))" for i in 1:D for j in i:D]
  str = join(terms, ", ")
  Meta.parse("SymTensorValue{$D,$S}(($str,))")
end

"""
    to_voigt(a::SymFourthOrderTensorValue{D,T,L}) -> TensorValue{M,M,T}

Encode a symmetric fourth-order tensor in Voigt notation as an M×M matrix,
where M = D(D+1)/2. Entry [I,J] of the result equals `a[i,j,k,l]`, where
(i,j) and (k,l) are the index pairs corresponding to Voigt indices I and J.
"""
@generated function to_voigt(a::SymFourthOrderTensorValue{D,T,L}) where {D,T,L}
  M = D*(D+1)÷2
  pairs = _voigt_pairs(D)
  terms = String[]
  for J in 1:M, I in 1:M
    i, j = pairs[I]; k, l = pairs[J]
    push!(terms, "a[$i,$j,$k,$l]")
  end
  str = join(terms, ", ")
  Meta.parse("TensorValue{$M,$M,T}(($str,))")
end

"""
    from_voigt(m::TensorValue{N,N,T,L}) -> SymFourthOrderTensorValue

Decode a Voigt matrix to a symmetric fourth-order tensor.
`N` must satisfy N = D(D+1)/2 for some integer D.
"""
@generated function from_voigt(m::TensorValue{N,N,T,L}) where {N,T,L}
  D = (isqrt(1+8*N)-1) ÷ 2
  @assert N == D*(D+1)÷2 "from_voigt: TensorValue size $N×$N is not a valid Voigt matrix size"
  inv_map = _voigt_inv(D)
  terms = String[]
  for i in 1:D, j in i:D, k in 1:D, l in k:D
    I = inv_map[i,j]; J = inv_map[k,l]
    push!(terms, "m[$I,$J]")
  end
  str = join(terms, ", ")
  Meta.parse("SymFourthOrderTensorValue{$D,T}(($str,))")
end

"""
    to_mandel(a::SymFourthOrderTensorValue{D,T,L}) -> TensorValue{M,M,S}

Encode a symmetric fourth-order tensor in Mandel notation as an M×M matrix,
where M = D(D+1)/2. Off-diagonal Voigt indices are scaled by √2 on each
axis, so double off-diagonal entries are scaled by 2 (exact).
The output eltype `S` is `promote_op(*, typeof(√2), T)`.
"""
@generated function to_mandel(a::SymFourthOrderTensorValue{D,T,L}) where {D,T,L}
  M = D*(D+1)÷2
  S = _mandel_eltype(T)
  pairs = _voigt_pairs(D)
  terms = String[]
  for J in 1:M, I in 1:M
    i, j = pairs[I]; k, l = pairs[J]
    if i==j && k==l
      push!(terms, "convert($S, a[$i,$j,$k,$l])")
    elseif i==j || k==l
      push!(terms, "sqrt($S(2)) * a[$i,$j,$k,$l]")
    else
      push!(terms, "$S(2) * a[$i,$j,$k,$l]")
    end
  end
  str = join(terms, ", ")
  Meta.parse("TensorValue{$M,$M,$S}(($str,))")
end

"""
    from_mandel(m::TensorValue{N,N,T,L}) -> SymFourthOrderTensorValue

Decode a Mandel matrix to a symmetric fourth-order tensor.
`N` must satisfy N = D(D+1)/2 for some integer D.
"""
@generated function from_mandel(m::TensorValue{N,N,T,L}) where {N,T,L}
  D = (isqrt(1+8*N)-1) ÷ 2
  @assert N == D*(D+1)÷2 "from_mandel: TensorValue size $N×$N is not a valid Mandel matrix size"
  S = _mandel_eltype(T)
  inv_map = _voigt_inv(D)
  terms = String[]
  for i in 1:D, j in i:D, k in 1:D, l in k:D
    I = inv_map[i,j]; J = inv_map[k,l]
    if i==j && k==l
      push!(terms, "convert($S, m[$I,$J])")
    elseif i==j || k==l
      push!(terms, "m[$I,$J] / sqrt($S(2))")
    else
      push!(terms, "m[$I,$J] / $S(2)")
    end
  end
  str = join(terms, ", ")
  Meta.parse("SymFourthOrderTensorValue{$D,$S}(($str,))")
end

