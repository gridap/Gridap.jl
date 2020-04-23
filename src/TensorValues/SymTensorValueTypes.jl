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

SymTensorValue(data::NTuple{L,T}) where {L,T}                = SymTensorValue{floor(Int,sqrt(L*2)),T}(data)
SymTensorValue{D}(data::NTuple{L,T}) where {D,L,T}           = SymTensorValue{D,T}(NTuple{L,T}(data))
SymTensorValue{D,T1}(data::NTuple{L,T2}) where {D,L,T1,T2}   = SymTensorValue{D,T1}(NTuple{L,T1}(data))
SymTensorValue{D,T1,L}(data::NTuple{L,T2}) where {D,L,T1,T2} = SymTensorValue{D,T1}(NTuple{L,T1}(data))

# SymTensorValue Vararg constructor

SymTensorValue(data::T...) where {T}              = SymTensorValue(Tuple(data))
SymTensorValue{D}(data::T...) where {D,T}         = SymTensorValue{D}(Tuple(data))
SymTensorValue{D,T1}(data::T2...) where {D,T1,T2} = SymTensorValue{D,T1}(Tuple(data))

# SymTensorValue single AbstractMatrix argument constructor

#From Square Matrices
_flatten_upper_triangle(data::AbstractArray,::Val{D}) where D = Tuple(data[i,j] for i in 1:D for j in i:D)

SymTensorValue(data::AbstractMatrix{T}) where {T} = ((D1,D2)=size(data); SymTensorValue{D1}(data))
SymTensorValue{D}(data::AbstractMatrix{T}) where {D,T} = SymTensorValue{D,T}(_flatten_upper_triangle(data,Val{D}()))
SymTensorValue{D,T1}(data::AbstractMatrix{T2}) where {D,T1,T2} = SymTensorValue{D,T1}(_flatten_upper_triangle(data,Val{D}()))
SymTensorValue{D,T1,L}(data::AbstractMatrix{T2}) where {D,T1,T2,L} = SymTensorValue{D,T1,L}(_flatten_upper_triangle(data,Val{D}()))

###############################################################
# Conversions (SymTensorValue)
###############################################################

function _SymTensorValue_to_array(arg::SymTensorValue{D,T,L}) where {D,T,L}
    z = zeros(MMatrix{D,D,T})
    for j in 1:D
        for i in j:D
            z[j,i] = z[i,j] = arg[i,j]
        end
    end
    z
end

# Direct conversion
convert(::Type{<:SymTensorValue{D,T}}, arg::AbstractArray) where {D,T} = SymTensorValue{D,T}(arg)
convert(::Type{<:SymTensorValue{D,T}}, arg::Tuple) where {D,T} = SymTensorValue{D,T}(arg)

# Inverse conversion
convert(::Type{<:MMatrix{D,D,T}}, arg::SymTensorValue) where {D,T} = _SymTensorValue_to_array(arg)
convert(::Type{<:SMatrix{D,D,T}}, arg::SymTensorValue) where {D,T} = SMatrix{D,D,T}(_SymTensorValue_to_array(arg))
convert(::Type{<:NTuple{L,T}}, arg::SymTensorValue) where {L,T} = NTuple{L,T}(Tuple(arg))

# Internal conversion
convert(::Type{<:SymTensorValue{D,T}}, arg::SymTensorValue{D}) where {D,T} = SymTensorValue{D,T}(Tuple(arg))
convert(::Type{<:SymTensorValue{D,T}}, arg::SymTensorValue{D,T}) where {D,T} = arg

###############################################################
# Other constructors and conversions (SymTensorValue)
###############################################################

zero(::Type{<:SymTensorValue{D,T}}) where {D,T} = (L=Int(D*(D+1)/2);SymTensorValue{D,T}(tfill(zero(T),Val{L}())))
zero(::Type{<:SymTensorValue{D,T,L}}) where {D,T,L} = SymTensorValue{D,T}(tfill(zero(T),Val{L}()))
zero(::SymTensorValue{D,T,L}) where {D,T,L} = zero(SymTensorValue{D,T,L})

@generated function one(::Type{<:SymTensorValue{D,T}}) where {D,T}
  str = join(["$i==$j ? one(T) : zero(T), " for i in 1:D for j in i:D])
  Meta.parse("SymTensorValue{D,T}(($str))")
end
one(::SymTensorValue{D,T}) where {D,T} = one(SymTensorValue{D,T})

mutable(::Type{<:SymTensorValue{D,T}}) where {D,T} = MMatrix{D,D,T}
mutable(::SymTensorValue{D,T}) where {D,T} = mutable(SymTensorValue{D,T})

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

n_components(::Type{<:SymTensorValue{D}}) where {D} = length(SymTensorValue{D})
n_components(::SymTensorValue{D}) where {D} = n_components(SymTensorValue{D})

