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

###############################################################
# Conversions (SymTensorValue)
###############################################################

function convert(TT::Type{<:Union{SymTensorValue,SymTensorValue{D,T1},SymTensorValue{D,T1,L}}}, 
                arg::
                    Union{
                        NTuple{L,T2},
                    }) where {D,T1,T2,L}
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

function convert(T::Type{<:Union{NTuple,NTuple{L,T1}}}, arg::SymTensorValue{D,T2,L}) where {D,T1,T2,L}
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
one(::SymTensorValue{D,T}) where {D,T} = one(TensorValue{D,T})

mutable(::Type{<:SymTensorValue{D,T}}) where {D,T} = MMatrix{D,D,T}
mutable(::SymTensorValue{D,T}) where {D,T} = mutable(SymTensorValue{D,T})

change_eltype(::Type{SymTensorValue{D,T1,L}},::Type{T2}) where {D,T1,T2,L} = SymTensorValue{D,T2,L}
change_eltype(::SymTensorValue{D,T1,L},::Type{T2}) where {D,T1,T2,L} = change_eltype(SymTensorValue{D,T1,L},T2)

function SymTensorValueToArray(arg::SymTensorValue{D,T,L}) where {D,T,L}
    z = zeros(T,D,D)
    vector = collect(Tuple(arg))
    for i in 1:D
        index = _getindex(arg,i)
        range = index:index+(D-i)
        z[i,i:D] = z[i:D,i] = vector[range]
    end
    z
end

SMatrix(arg::SymTensorValue{D,T,L}) where {D,T,L} = SMatrix{D,D,T,L}(SymTensorValueToArray(arg.data))
SArray(arg::SymTensorValue{D,T,L}) where {D,T,L} =  SMatrix(arg)
get_array(arg::SymTensorValue{D,T,L}) where {D,T,L} = SMatrix(arg)

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

