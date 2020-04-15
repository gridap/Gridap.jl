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

function _getindex(arg::SymTensorValue{D},i::Integer) where {D}
    index=((i-1)*D)-sum(1:i-1)+i
end

function _getindex(arg::SymTensorValue{D},i::Integer,j::Integer) where {D}
    _j,_i=sort([i,j])
    index=((_j-1)*D)-sum(1:_j-1)+_i
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
SymTensorValue{D}(data::NTuple{L,T}) where {D,L,T}           = SymTensorValue{D,T}(data)
SymTensorValue{D,T1}(data::NTuple{L,T2}) where {D,L,T1,T2}   = SymTensorValue{D,T1}(NTuple{L,T1}(data))
SymTensorValue{D,T1,L}(data::NTuple{L,T2}) where {D,L,T1,T2} = SymTensorValue{D,T1}(NTuple{L,T1}(data))

# SymTensorValue Vararg constructor

SymTensorValue(data::Real...)                  = (L=length(data);SymTensorValue{floor(Int,sqrt(L*2))}(NTuple{L}(data)))
SymTensorValue{D}(data::Real...) where {D}     = (L=length(data);SymTensorValue{D}(NTuple{L}(data)))
SymTensorValue{D,T}(data::Real...) where {D,T} = (L=length(data);SymTensorValue{D,T}(NTuple{L,T}(data)))

#From Square Matrices
# SymTensorValue -> b[[(D*(i-1))+j for i in 1:D for j in i:D]]

# SymTensorValue single SVector, MVector, SMatrix, MMatrix and AbstractMatrix argument constructor

function _FlattenUpperTriangle(data::
                Union{
                    SMatrix{D1,D2,T2,L},
                    MMatrix{D1,D2,T2,L},
                    AbstractMatrix{T2}
                }) where {D1,D2,T1,T2,L}
    PD1 = (@isdefined D1) ? D1 : size(data)[1]
    PD2 = (@isdefined D2) ? D2 : size(data)[2]
    [data[i,j] for i in 1:PD1 for j in i:PD2]
end

function SymTensorValue(data::
                Union{
                    SMatrix{D,D,T2,L},
                    MMatrix{D,D,T2,L},
                    AbstractMatrix{T2}
                }) where {D,T1,T2,L}
    PD1 = (@isdefined D) ? D : size(data)[1]
    PD2 = (@isdefined D) ? D : size(data)[2]
    @assert PD1 == PD2
    ut=_FlattenUpperTriangle(data)
    SymTensorValue{PD1,T2}(NTuple{length(ut),T2}(ut))
end

function SymTensorValue{D}(data::
                Union{
                    SMatrix{D,D,T2,L},
                    MMatrix{D,D,T2,L},
                    AbstractMatrix{T2}
                }) where {D,T1,T2,L}
    ut=_FlattenUpperTriangle(data)
    SymTensorValue{D,T2}(NTuple{length(ut),T2}(ut))
end

function SymTensorValue{D,T1}(data::
                Union{
                    SMatrix{D,D,T2,L},
                    MMatrix{D,D,T2,L},
                    AbstractMatrix{T2}
                }) where {D,T1,T2,L}
    ut=_FlattenUpperTriangle(data)
    SymTensorValue{D,T1}(NTuple{length(ut),T1}(ut))
end

function SymTensorValue{D,T1,L1}(data::
                Union{
                    SMatrix{D,D,T2,L2},
                    MMatrix{D,D,T2,L2},
                    AbstractMatrix{T2}
                }) where {D,T1,T2,L1,L2}
    ut=_FlattenUpperTriangle(data)
    SymTensorValue{D,T1}(NTuple{L1,T1}(ut))
end

###############################################################
# Conversions (SymTensorValue)
###############################################################

function SymTensorValueToArray(arg::SymTensorValue{D,T,L}) where {D,T,L}
    z = zeros(T,D,D)
    data = collect(Tuple(arg))
    for i in 1:D
        index = _getindex(arg,i)
        range = index:index+(D-i)
        z[i,i:D] = z[i:D,i] = data[range]
    end
    z
end

function convert(::Type{<:Union{SymTensorValue,SymTensorValue{D,T1},SymTensorValue{D,T1,L}}}, 
                arg::NTuple{L,T2}) where {D,T1,T2,L}
    PT = (@isdefined T1) ? T1 : T2
    SymTensorValue{D,PT}(arg)
end

function convert(::Type{Union{SMatrix,SMatrix{D,D}}}, arg::SymTensorValue{D,T2,L2}) where {D,T1,T2,L1,L2}
    PT = (@isdefined T1) ? T1 : T2
    SMatrix{D,D,PT}(SymTensorValueToArray(arg))
end

function convert(::Type{<:SMatrix{D,D,T1}}, arg::SymTensorValue{D,T2,L}) where {D,T1,T2,L}
    SMatrix{D,D,T1}(SymTensorValueToArray(arg))
end

function convert(::Type{<:SMatrix{D,D}}, arg::SymTensorValue{D,T2,L}) where {D,T2,L}
    SMatrix{D,D,T2}(SymTensorValueToArray(arg))
end

function convert(::Type{Union{MMatrix,MMatrix{D,D,T1,L1}}}, arg::SymTensorValue{D,T2,L2}) where {D,T1,T2,L1,L2}
    PT = (@isdefined T1) ? T1 : T2
    MMatrix{D,D,PT}(SymTensorValueToArray(arg))
end

function convert(::Type{<:MMatrix{D,D,T1}}, arg::SymTensorValue{D,T2,L}) where {D,T1,T2,L}
    MMatrix{D,D,T1}(SymTensorValueToArray(arg))
end

function convert(::Type{<:MMatrix{D,D}}, arg::SymTensorValue{D,T2,L}) where {D,T2,L}
    MMatrix{D,D,T2}(SymTensorValueToArray(arg))
end

function convert(::Type{<:Union{NTuple,NTuple{L,T1}}}, arg::SymTensorValue{D,T2,L}) where {D,T1,T2,L}
    PT = (@isdefined T1) ? T1 : T2
    NTuple{L,PT}(arg.data)
end

###############################################################
# Other constructors and conversions (SymTensorValue)
###############################################################

zero(::Type{<:SymTensorValue{D,T}}) where {D,T} = (L=Int(D*(D+1)/2);SymTensorValue{D,T}(NTuple{L,T}(zeros(T,L))))
zero(::SymTensorValue{D,T}) where {D,T} = zero(SymTensorValue{D,T})

@generated function one(::Type{<:SymTensorValue{D,T}}) where {D,T}
  str = join(["$i==$j ? one(T) : zero(T), " for i in 1:D for j in i:D])
  Meta.parse("SymTensorValue{D,T}(($str))")
end
one(::SymTensorValue{D,T}) where {D,T} = one(SymTensorValue{D,T})

mutable(::Type{<:SymTensorValue{D,T}}) where {D,T} = MMatrix{D,D,T}
mutable(::SymTensorValue{D,T}) where {D,T} = mutable(SymTensorValue{D,T})

change_eltype(::Type{SymTensorValue{D,T1,L}},::Type{T2}) where {D,T1,T2,L} = SymTensorValue{D,T2,L}
change_eltype(::SymTensorValue{D,T1,L},::Type{T2}) where {D,T1,T2,L} = change_eltype(SymTensorValue{D,T1,L},T2)

SMatrix(arg::SymTensorValue{D,T,L}) where {D,T,L} = convert(SMatrix{D,D,T}, arg)
SArray(arg::SymTensorValue{D,T,L}) where {D,T,L} =  convert(SMatrix{D,D,T}, arg)
get_array(arg::SymTensorValue{D,T,L}) where {D,T,L} = convert(SMatrix{D,D,T}, arg)

###############################################################
# Introspection (SymTensorValue)
###############################################################


eltype(::Type{<:SymTensorValue{D,T}}) where {D,T} = T
eltype(::SymTensorValue{D,T}) where {D,T} = eltype(SymTensorValue{D,T})

size(::Type{SymTensorValue{D}}) where {D} = (D,D)
size(::Type{SymTensorValue{D,T}}) where {D,T} = (D,D)
size(::Type{SymTensorValue{D,T,L}}) where {D,T,L} = (D,D)
size(::SymTensorValue{D,T}) where {D,T} = size(SymTensorValue{D,T})

length(::Type{SymTensorValue{D}}) where {D} = Int(D*(D+1)/2)
length(::Type{SymTensorValue{D,T}}) where {D,T} = length(SymTensorValue{D})
length(::Type{SymTensorValue{D,T,L}}) where {D,T,L} = L
length(::SymTensorValue{D,T,L}) where {D,T,L} = length(SymTensorValue{D,T,L})

n_components(::Type{SymTensorValue{D}}) where {D} = length(SymTensorValue{D})
n_components(::Type{SymTensorValue{D,T,L}}) where {D,T,L} = length(SymTensorValue{D,T,L})
n_components(::SymTensorValue{D,T,L}) where {D,T,L} = n_components(SymTensorValue{D,T,L})

