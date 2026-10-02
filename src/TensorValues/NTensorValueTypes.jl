
"""
    struct NTensorValue{S,T,N,L} <: MultiValue{S,T,N,L}

Type representing a `N`th tensor of size `S` with no dependant components.
Thus, it always hold `L` = `prod(S)` components.

If `Val(N)` is given as first argument of the constructor, `allequal(S)` is assumed.

`NTensorValue` is a `MultiValue` wrapper for `SArray{S,T,N,L}` with `N` ≥ 3.
"""
struct NTensorValue{S,T,N,L} <: MultiValue{S,T,N,L} #<:NTuple{N}
  data::NTuple{L,T}

  function NTensorValue(data::SArray{S,T,N,L}) where {S,T,N,L}
    @check N > 2 "NTensorValue starts from 3th order, use VectorValue or TensorValue instead"
    new{S,T,N,L}(NTuple{L,T}(data))
  end

  function NTensorValue{S,T,N}(data::NTuple{L,T}) where {S,T,N,L}
    @check N > 2 "NTensorValue starts from 3th order, use VectorValue or TensorValue instead"
    new{S,T,N,L}(data)
  end
end

@inline Base.Tuple(arg::NTensorValue) = Tuple(arg.data)

function promote_rule(
  ::Type{<:NTensorValue{S,Ta,N}},
  ::Type{<:NTensorValue{S,Tb,N}}) where {N,S,Ta,Tb}

  T = promote_type(Ta,Tb)
  NTensorValue{S,T,N}
end

# Empty NTensorValue constructor

NTensorValue(::Val{N}) where N       = NTensorValue{NTuple{N,0},Int,N}(NTuple{0,Int}())
NTensorValue(::Val{N}, data::NTuple{0}) where N = NTensorValue(Val(N))
NTensorValue{S}(data::NTuple{0}) where S = NTensorValue{S,Int}(data)

# NTensorValue single NTuple argument constructor

@generated function NTensorValue(::Val{N}, data::NTuple{L,T}) where {N,L,T}
  msg = "Invalid number of scalar arguments in NTensorValue{$N} constructor, expecting a $(N)th power number, got L=$L"
  D = round(Int, L^(1/N))
  @assert L == D^N msg
  S = NTuple{N,D}
  quote
    NTensorValue{$S,T,N}(data)
  end
end
NTensorValue(::Val{N}, data...) where N = NTensorValue(Val(N),data)
NTensorValue(::Val{N}, data::AbstractArray{T,N}) where {N,T} = NTensorValue(Val(N), data...)

@generated NTensorValue{S}(data::NTuple{L,T}) where {S,L,T}   = :( NTensorValue{S,T,$(length(S.parameters)),}(data) )
@generated NTensorValue{S,T}(data::NTuple{L,T}) where {S,L,T} = :( NTensorValue{S,T,$(length(S.parameters)),}(data) )

# NTensorValue single Tuple argument constructor

#NTensorValue(data::Tuple) = NTensorValue(promote(data...))
#NTensorValue{D}(data::Tuple) where {D} = NTensorValue{D,D,D}(promote(data...))
NTensorValue(::Val{N}, data::Tuple) where N = NTensorValue(Val(N), promote(data...))
NTensorValue{S}(data::Tuple) where S = NTensorValue{S}(promote(data...))
NTensorValue{S,T1}(data::Tuple) where {S,T1} = NTensorValue{S,T1}(NTuple{length(data),T1}(data))
NTensorValue{S,T1,N}(data::Tuple) where {S,T1,N} = NTensorValue{S,T1,N}(NTuple{length(data),T1}(data))
NTensorValue{S,T1,N,L}(data::Tuple) where {S,T1,N,L} = NTensorValue{S,T1,N}(NTuple{L,T1}(data))

# NTensorValue Vararg constructor

#NTensorValue(data::Number...) = NTensorValue(data)
NTensorValue{S}(data::Number...) where S =  NTensorValue{S}(data)
NTensorValue{S,T}(data::Number...) where {S,T} =  NTensorValue{S,T}(data)
NTensorValue{S,T,N}(data::Number...) where {S,T,N} =  NTensorValue{S,T,N}(data)
NTensorValue{S,T,N,L}(data::Number...) where {S,T,N,L} =  NTensorValue{S,T,N,L}(data)

# NTensorValue single AbstractArray{3,T} argument constructor

function NTensorValue(data::AbstractArray{T,N}) where {T,N}
  s = size(data)
  S = Tuple{s...}
  L = length(data)
  NTensorValue{S,T,N}(NTuple{L,T}(data))
end
NTensorValue{S}(data::AbstractArray{T,N}) where {S,T,N} = (L=length(data);NTensorValue{S,T,N}(NTuple{L,T}(data)))
NTensorValue{S,T1}(data::AbstractArray{T2,N}) where {S,T1,T2,N} = (L=length(data);NTensorValue{S,T1,N}(NTuple{L,T1}(data)))
NTensorValue{S,T1,N}(data::AbstractArray{T2,N}) where {S,T1,T2,N} = (L=length(data);NTensorValue{S,T1,N}(NTuple{L,T1}(data)))
NTensorValue{S,T1,N,L}(data::AbstractArray{T2,N}) where {S,T1,T2,N,L} = NTensorValue{S,T1,N}(NTuple{L,T1}(data))

###############################################################
# Conversions (NTensorValue)
###############################################################

# Inverse conversion
convert(::Type{<:SArray{S,T}}, arg::NTensorValue) where {S,T} = SArray{S,T}(Tuple(arg))
convert(::Type{<:MArray{S,T}}, arg::NTensorValue) where {S,T} = MArray{S,T}(Tuple(arg))

# Internal conversion
convert(::Type{<:NTensorValue{S,T}}, arg::NTensorValue{S}) where {S,T} = NTensorValue{S,T}(Tuple(arg))
convert(::Type{<:NTensorValue{S,T,N}}, arg::NTensorValue{S,T,N}) where {S,T,N} = arg

# Construction from SArray or MArray
function MultiValue(a::StaticArray{S,T,N}) where {T,N,S<:Tuple}
  # All 0-2 order tensors must implement this method explicitly
  @check N > 2 "A specialization of `MultiValue(::StaticArrays)` is missing for $(typeof(a))"
  convert(NTensorValue{S,T,N}, a)
end

# other
#change_eltype(::Type{<:NTensorValue{D1,D2,D3}},::Type{T2}) where {D1,D2,D3,D4,T2} = NTensorValue{D1,D2,D3,D4,T2}
change_eltype(::Type{NTensorValue{S,T1,N,L}},::Type{T2}) where {S,T1,T2,N,L} = NTensorValue{S,T2,N,L}

