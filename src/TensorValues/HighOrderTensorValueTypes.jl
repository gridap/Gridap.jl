
"""
    struct HighOrderTensorValue{S,T,N,L} <: MultiValue{S,T,N,L}

Type representing a `N`th tensor of size `S` with no dependant components.
Thus, it always hold `L` = `prod(S)` components.

If `Val(N)` is given as first argument of the constructor, `allequal(S)` is assumed.

`HighOrderTensorValue` is a `MultiValue` wrapper for `SArray{S,T,N,L}` with `N` ≥ 3.
"""
struct HighOrderTensorValue{S,T,N,L} <: MultiValue{S,T,N,L} #<:NTuple{N}
  data::SArray{S,T,N,L}

  function HighOrderTensorValue(data::SArray{S,T,N,L}) where {S,T,N,L}
    @check N > 2 "HighOrderTensorValue starts from 3th order, use VectorValue or TensorValue instead"
    new{S,T,N,L}(data)
  end

  function HighOrderTensorValue{S,T,N}(data::NTuple{L,T}) where {S,T,N,L}
    @check N > 2 "HighOrderTensorValue starts from 3th order, use VectorValue or TensorValue instead"
    new{S,T,N,L}(data)
  end
end

@inline Base.Tuple(arg::HighOrderTensorValue) = Tuple(arg.data)

function promote_rule(
  ::Type{<:HighOrderTensorValue{S,Ta,N}},
  ::Type{<:HighOrderTensorValue{S,Tb,N}}) where {N,S,Ta,Tb}

  T = promote_type(Ta,Tb)
  HighOrderTensorValue{S,T,N}
end

# Empty HighOrderTensorValue constructor

HighOrderTensorValue()                       = HighOrderTensorValue{4,Tuple{0,0,0,0},Int}(NTuple{0,Int}())
HighOrderTensorValue(::Val{N}) where N       = HighOrderTensorValue{Tuple{ntuple(_->0,Val(N))...},Int,N}(NTuple{0,Int}())
HighOrderTensorValue(::Val{N}, data::NTuple{0}) where N = HighOrderTensorValue(Val(N))
HighOrderTensorValue{S}(data::NTuple{0}) where S = HighOrderTensorValue{S,Int}(data)

# HighOrderTensorValue single NTuple argument constructor

#HighOrderTensorValue(::Val{N}, data::NTuple{L,T}) where {N,L,T} = HighOrderTensorValue{N}(data)
#@generated function HighOrderTensorValue{N}(data::NTuple{L,T}) where {N,L,T}
@generated function HighOrderTensorValue(::Val{N}, data::NTuple{L,T}) where {N,L,T}
  msg = "Invalid number of scalar arguments in HighOrderTensorValue{$N} constructor, expecting a $(N)th power number, got L=$L"
  V = L^(1/N)
  @assert floor(Int,V) == ceil(Int,V) msg
  D=Int(V)
  S = Tuple{ntuple(_->D, Val(N))...}
  quote
    HighOrderTensorValue{$S,T,N}(data)
  end
end
HighOrderTensorValue(::Val{N}, data...) where N = HighOrderTensorValue(Val(N),data)
HighOrderTensorValue(::Val{N}, data::AbstractArray{T,N}) where {N,T} = HighOrderTensorValue(Val(N), data...)

@generated HighOrderTensorValue{S}(data::NTuple{L,T}) where {S,L,T}   = :( HighOrderTensorValue{S,T,$(length(S.parameters)),}(data) )
@generated HighOrderTensorValue{S,T}(data::NTuple{L,T}) where {S,L,T} = :( HighOrderTensorValue{S,T,$(length(S.parameters)),}(data) )

# HighOrderTensorValue single Tuple argument constructor

#HighOrderTensorValue(data::Tuple) = HighOrderTensorValue(promote(data...))
#HighOrderTensorValue{D}(data::Tuple) where {D} = HighOrderTensorValue{D,D,D}(promote(data...))
HighOrderTensorValue(::Val{N}, data::Tuple) where N = HighOrderTensorValue(Val(N), promote(data...))
HighOrderTensorValue{S}(data::Tuple) where S = HighOrderTensorValue{S}(promote(data...))
HighOrderTensorValue{S,T1}(data::Tuple) where {S,T1} = HighOrderTensorValue{S,T1}(NTuple{length(data),T1}(data))
HighOrderTensorValue{S,T1,N}(data::Tuple) where {S,T1,N} = HighOrderTensorValue{S,T1,N}(NTuple{length(data),T1}(data))
HighOrderTensorValue{S,T1,N,L}(data::Tuple) where {S,T1,N,L} = HighOrderTensorValue{S,T1,N}(NTuple{L,T1}(data))

# HighOrderTensorValue Vararg constructor

#HighOrderTensorValue(data::Number...) = HighOrderTensorValue(data)
HighOrderTensorValue{S}(data::Number...) where S =  HighOrderTensorValue{S}(data)
HighOrderTensorValue{S,T}(data::Number...) where {S,T} =  HighOrderTensorValue{S,T}(data)
HighOrderTensorValue{S,T,N}(data::Number...) where {S,T,N} =  HighOrderTensorValue{S,T,N}(data)
HighOrderTensorValue{S,T,N,L}(data::Number...) where {S,T,N,L} =  HighOrderTensorValue{S,T,N,L}(data)

# HighOrderTensorValue single AbstractArray{3,T} argument constructor

function HighOrderTensorValue(data::AbstractArray{T,N}) where {T,N}
  s = size(data)
  S = Tuple{s...}
  L = length(data)
  HighOrderTensorValue{S,T,N}(NTuple{L,T}(data))
end
HighOrderTensorValue{S}(data::AbstractArray{T,N}) where {S,T,N} = (L=length(data);HighOrderTensorValue{S,T,N}(NTuple{L,T}(data)))
HighOrderTensorValue{S,T1}(data::AbstractArray{T2,N}) where {S,T1,T2,N} = (L=length(data);HighOrderTensorValue{S,T1,N}(NTuple{L,T1}(data)))
HighOrderTensorValue{S,T1,N}(data::AbstractArray{T2,N}) where {S,T1,T2,N} = (L=length(data);HighOrderTensorValue{S,T1,N}(NTuple{L,T1}(data)))
HighOrderTensorValue{S,T1,N,L}(data::AbstractArray{T2,N}) where {S,T1,T2,N,L} = HighOrderTensorValue{S,T1,N}(NTuple{L,T1}(data))

###############################################################
# Conversions (HighOrderTensorValue)
###############################################################

# Inverse conversion
convert(::Type{<:SArray{S,T}}, arg::HighOrderTensorValue) where {S,T} = SArray{S,T}(Tuple(arg))
convert(::Type{<:MArray{S,T}}, arg::HighOrderTensorValue) where {S,T} = MArray{S,T}(Tuple(arg))

# Internal conversion
convert(::Type{<:HighOrderTensorValue{S,T}}, arg::HighOrderTensorValue{S}) where {S,T} = HighOrderTensorValue{S,T}(Tuple(arg))
convert(::Type{<:HighOrderTensorValue{S,T,N}}, arg::HighOrderTensorValue{S,T,N}) where {S,T,N} = arg

# Construction from SArray or MArray
function MultiValue(a::StaticArray{S,T,N}) where {T,N,S<:Tuple}
  # All 0-2 order tensors must implement this method explicitly
  @check N > 2 "A specialization of `MultiValue(::StaticArrays)` is missing for $(typeof(a))"
  convert(HighOrderTensorValue{S,T,N}, a)
end

# other
#change_eltype(::Type{<:HighOrderTensorValue{D1,D2,D3}},::Type{T2}) where {D1,D2,D3,D4,T2} = HighOrderTensorValue{D1,D2,D3,D4,T2}
change_eltype(::Type{HighOrderTensorValue{S,T1,N,L}},::Type{T2}) where {S,T1,T2,N,L} = HighOrderTensorValue{S,T2,N,L}

