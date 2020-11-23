###############################################################
# SymTensorValue Type
###############################################################

"""
Type representing a symmetric second-order tensor
"""
struct SymTensorValue{D,T,L} <: MultiValue{Tuple{D,D},T,2,L}
    data::NTuple{L,T}
    function SymTensorValue{D,T}(data::NTuple{L,T}) where {D,T,L}
        @assert L == D*(D+1)/2
        new{D,T,L}(data)
    end
end

###############################################################
# Constructors (SymTensorValue)
###############################################################

# Empty SymTensorValue constructor

SymTensorValue()                   = SymTensorValue{0,Int}(NTuple{0,Int}())
SymTensorValue{0}() where {T}      = SymTensorValue{0,Int}(NTuple{0,Int}())
SymTensorValue{0,T}() where {T}    = SymTensorValue{0,T}(NTuple{0,T}())
SymTensorValue(data::NTuple{0})    = SymTensorValue{0,Int}(data)
SymTensorValue{0}(data::NTuple{0}) = SymTensorValue{0,Int}(data)

# SymTensorValue single NTuple argument constructor

@generated function SymTensorValue(data::NTuple{L,T}) where {L,T}
  msg = "Invalid number of scalar arguments in SymTensorValue constructor"
  V = (sqrt(1+8*L)-1)/2
  @assert floor(Int,V) == ceil(Int,V) msg
  D = Int(V)
  quote
    SymTensorValue{$D,T}(data)
  end
end
SymTensorValue{D}(data::NTuple{L,T}) where {D,L,T}           = SymTensorValue{D,T}(data)
SymTensorValue{D,T1}(data::NTuple{L,T2}) where {D,L,T1,T2}   = SymTensorValue{D,T1}(NTuple{L,T1}(data))
SymTensorValue{D,T1,L}(data::NTuple{L,T2}) where {D,L,T1,T2} = SymTensorValue{D,T1}(NTuple{L,T1}(data))

# SymTensorValue single Tuple argument constructor

SymTensorValue(data::Tuple)                      = SymTensorValue(promote(data...))
SymTensorValue{D}(data::Tuple) where {D}         = SymTensorValue{D}(promote(data...))
SymTensorValue{D,T1}(data::Tuple) where {D,T1}   = SymTensorValue{D,T1}(NTuple{length(data),T1}(data))
SymTensorValue{D,T1,L}(data::Tuple) where {D,T1,L}   = SymTensorValue{D,T1}(NTuple{L,T1}(data))

# SymTensorValue Vararg constructor

SymTensorValue(data...)                    = SymTensorValue(data)
SymTensorValue{D}(data...)    where {D}    = SymTensorValue{D}(data)
SymTensorValue{D,T1}(data...) where {D,T1} = SymTensorValue{D,T1}(data)
SymTensorValue{D,T1,L}(data...) where {D,T1,L} = SymTensorValue{D,T1}(data)

# SymTensorValue single AbstractMatrix argument constructor

#From Square Matrices
@generated function _flatten_upper_triangle(data::AbstractArray,::Val{D}) where D
  str = ""
  for i in 1:D
    for j in i:D
      str *= "data[i,j], "
    end
  end
  Meta.parse("($str)")
end

SymTensorValue(data::AbstractMatrix{T}) where {T} = ((D1,D2)=size(data); SymTensorValue{D1}(data))
SymTensorValue{D}(data::AbstractMatrix{T}) where {D,T} = SymTensorValue{D,T}(_flatten_upper_triangle(data,Val{D}()))
SymTensorValue{D,T1}(data::AbstractMatrix{T2}) where {D,T1,T2} = SymTensorValue{D,T1}(_flatten_upper_triangle(data,Val{D}()))
SymTensorValue{D,T1,L}(data::AbstractMatrix{T2}) where {D,T1,T2,L} = SymTensorValue{D,T1,L}(_flatten_upper_triangle(data,Val{D}()))

###############################################################
# Conversions (SymTensorValue)
###############################################################

@generated function _SymTensorValue_to_array(arg::SymTensorValue{D,T,L}) where {D,T,L}
  str = ""
  for j in 1:D
    for i in 1:D
      p = _2d_sym_tensor_linear_index(D,i,j)
      str *= "arg.data[$p], "
    end
  end
  Meta.parse("SMatrix{D,D,T}(($str))")
end

# Direct conversion
convert(::Type{<:SymTensorValue{D,T}}, arg::AbstractArray) where {D,T} = SymTensorValue{D,T}(arg)
convert(::Type{<:SymTensorValue{D,T}}, arg::Tuple) where {D,T} = SymTensorValue{D,T}(arg)

# Inverse conversion
convert(::Type{<:MMatrix{D,D,T}}, arg::SymTensorValue) where {D,T} = MMatrix{D,D,T}(_SymTensorValue_to_array(arg))
convert(::Type{<:SMatrix{D,D,T}}, arg::SymTensorValue) where {D,T} = _SymTensorValue_to_array(arg)
convert(::Type{<:NTuple{L,T}}, arg::SymTensorValue) where {L,T} = NTuple{L,T}(Tuple(arg))

# Internal conversion
convert(::Type{<:SymTensorValue{D,T}}, arg::SymTensorValue{D}) where {D,T} = SymTensorValue{D,T}(Tuple(arg))
convert(::Type{<:SymTensorValue{D,T}}, arg::SymTensorValue{D,T}) where {D,T} = arg

###############################################################
# Other constructors and conversions (SymTensorValue)
###############################################################

@generated function zero(::Type{<:SymTensorValue{D,T}}) where {D,T}
  L=Int(D*(D+1)/2)
  quote
    SymTensorValue{D,T}(tfill(zero(T),Val{$L}()))
  end
end

zero(::Type{<:SymTensorValue{D,T,L}}) where {D,T,L} = SymTensorValue{D,T}(tfill(zero(T),Val{L}()))
zero(::SymTensorValue{D,T,L}) where {D,T,L} = zero(SymTensorValue{D,T,L})

@generated function one(::Type{<:SymTensorValue{D,T}}) where {D,T}
  str = join(["$i==$j ? one(T) : zero(T), " for i in 1:D for j in i:D])
  Meta.parse("SymTensorValue{D,T}(($str))")
end
one(::SymTensorValue{D,T}) where {D,T} = one(SymTensorValue{D,T})

Mutable(::Type{<:SymTensorValue{D,T}}) where {D,T} = MMatrix{D,D,T}
Mutable(::SymTensorValue{D,T}) where {D,T} = Mutable(SymTensorValue{D,T})
mutable(a::SymTensorValue{D}) where D = MMatrix{D,D}(Tuple(get_array(a)))

change_eltype(::Type{SymTensorValue{D,T1,L}},::Type{T2}) where {D,T1,T2,L} = SymTensorValue{D,T2,L}
change_eltype(::SymTensorValue{D,T1,L},::Type{T2}) where {D,T1,T2,L} = change_eltype(SymTensorValue{D,T1,L},T2)

get_array(arg::SymTensorValue{D,T,L}) where {D,T,L} = convert(SMatrix{D,D,T}, arg)

###############################################################
# Introspection (SymTensorValue)
###############################################################

eltype(::Type{<:SymTensorValue{D,T}}) where {D,T} = T
eltype(::SymTensorValue{D,T}) where {D,T} = eltype(SymTensorValue{D,T})

size(::Type{<:SymTensorValue{D}}) where {D} = (D,D)
size(::SymTensorValue{D}) where {D} = size(SymTensorValue{D})

length(::Type{<:SymTensorValue{D}}) where {D} = D*D
length(::SymTensorValue{D}) where {D} = length(SymTensorValue{D})

num_components(::Type{<:SymTensorValue{D}}) where {D} = length(SymTensorValue{D})
num_components(::SymTensorValue{D}) where {D} = num_components(SymTensorValue{D})

