###############################################################
# SymTensorValue Type
###############################################################

"""
Type representing a symmetric fourth-order tensor
"""
struct SymFourthOrderTensorValue{D,T,L} <: MultiValue{Tuple{D,D,D,D},T,4,L}
  data::NTuple{L,T}
  function SymFourthOrderTensorValue{D,T}(data::NTuple{L,T}) where {D,T,L}
    @assert L == (D*(D+1)/2)^2
    new{D,T,L}(data)
  end
end

#function _getindex(arg::SymTensorValue{D},i::Integer) where {D}
#    index=((i-1)*D)-sum(1:i-1)+i
#end

#function _getindex(arg::SymTensorValue{D},i::Integer,j::Integer) where {D}
#    _j,_i=sort([i,j])
#    index=((_j-1)*D)-sum(1:_j-1)+_i
#end



###############################################################
# Constructors (SymTensorValue)
###############################################################

# Empty SymTensorValue constructor

SymFourthOrderTensorValue()                   = SymFourthOrderTensorValue{0,Int}(NTuple{0,Int}())
SymFourthOrderTensorValue{0}() where {T}      = SymFourthOrderTensorValue{0,Int}(NTuple{0,Int}())
SymFourthOrderTensorValue{0,T}() where {T}    = SymFourthOrderTensorValue{0,T}(NTuple{0,T}())
SymFourthOrderTensorValue(data::NTuple{0})    = SymFourthOrderTensorValue{0,Int}(data)
SymFourthOrderTensorValue{0}(data::NTuple{0}) = SymFourthOrderTensorValue{0,Int}(data)

# SymTensorValue single NTuple argument constructor

SymFourthOrderTensorValue(data::NTuple{L,T}) where {L,T}                = SymFourthOrderTensorValue{floor(Int,sqrt(sqrt(L*2))),T}(data)
SymFourthOrderTensorValue{D}(data::NTuple{L,T}) where {D,L,T}           = SymFourthOrderTensorValue{D,T}(data)
SymFourthOrderTensorValue{D,T1}(data::NTuple{L,T2}) where {D,L,T1,T2}   = SymFourthOrderTensorValue{D,T1}(NTuple{L,T1}(data))
SymFourthOrderTensorValue{D,T1,L}(data::NTuple{L,T2}) where {D,L,T1,T2} = SymFourthOrderTensorValue{D,T1}(NTuple{L,T1}(data))

# SymTensorValue Vararg constructor

SymFourthOrderTensorValue(data::Real...)                  = (L=length(data);SymFourthOrderTensorValue{floor(Int,sqrt(sqrt(L*2)))}(NTuple{L}(data)))
SymFourthOrderTensorValue{D}(data::Real...) where {D}     = (L=length(data);SymFourthOrderTensorValue{D}(NTuple{L}(data)))
SymFourthOrderTensorValue{D,T}(data::Real...) where {D,T} = (L=length(data);SymFourthOrderTensorValue{D,T}(NTuple{L,T}(data)))

###############################################################
# Conversions (SymTensorValue)
###############################################################

function convert(::Type{<:Union{SymFourthOrderTensorValue,
                                SymFourthOrderTensorValue{D,T1},
                                SymFourthOrderTensorValue{D,T1,L}}}, 
                arg::NTuple{L,T2}) where {D,T1,T2,L}
    PT = (@isdefined T1) ? T1 : T2
    SymFourthOrderTensorValue{D,PT}(arg)
end

function convert(::Type{<:Union{NTuple,NTuple{L,T1}}}, 
                arg::SymFourthOrderTensorValue{D,T2,L}) where {D,T1,T2,L}
    PT = (@isdefined T1) ? T1 : T2
    NTuple{L,PT}(arg.data)
end

###############################################################
# Other constructors and conversions (SymTensorValue)
###############################################################

zero(::Type{<:SymFourthOrderTensorValue{D,T}}) where {D,T} = (L=Int((D*(D+1)/2)^2);SymFourthOrderTensorValue{D,T}(NTuple{L,T}(zeros(T,L))))
zero(::SymFourthOrderTensorValue{D,T}) where {D,T} = zero(SymFourthOrderTensorValue{D,T})

#@generated function one(::Type{<:SymFourthOrderTensorValue{D,T}}) where {D,T}
#  str = join(["$i==$j ? one(T) : zero(T), " for i in 1:D for j in i:D])
#  Meta.parse("SymFourthOrderTensorValue{D,T}(($str))")
#end
#one(::SymFourthOrderTensorValue{D,T}) where {D,T} = one(TensorValue{D,T})

change_eltype(::Type{SymFourthOrderTensorValue{D,T1,L}},::Type{T2}) where {D,T1,T2,L} = SymFourthOrderTensorValue{D,T2,L}
change_eltype(::SymFourthOrderTensorValue{D,T1,L},::Type{T2}) where {D,T1,T2,L} = change_eltype(SymFourthOrderTensorValue{D,T1,L},T2)

###############################################################
# Introspection (SymTensorValue)
###############################################################


eltype(::Type{<:SymFourthOrderTensorValue{D,T}}) where {D,T} = T
eltype(::SymFourthOrderTensorValue{D,T}) where {D,T} = eltype(SymFourthOrderTensorValue{D,T})

size(::Type{SymFourthOrderTensorValue{D}}) where {D} = (D,D,D,D)
size(::Type{SymFourthOrderTensorValue{D,T}}) where {D,T} = (D,D,D,D)
size(::Type{SymFourthOrderTensorValue{D,T,L}}) where {D,T,L} = (D,D,D,D)
size(::SymFourthOrderTensorValue{D,T}) where {D,T} = size(SymFourthOrderTensorValue{D,T})

length(::Type{SymFourthOrderTensorValue{D}}) where {D} = Int((D*(D+1)/2)^2)
length(::Type{SymFourthOrderTensorValue{D,T}}) where {D,T} = length(SymFourthOrderTensorValue{D})
length(::Type{SymFourthOrderTensorValue{D,T,L}}) where {D,T,L} = L
length(::SymFourthOrderTensorValue{D,T,L}) where {D,T,L} = length(SymFourthOrderTensorValue{D,T,L})

n_components(::Type{SymFourthOrderTensorValue{D}}) where {D} = length(SymFourthOrderTensorValue{D})
n_components(::Type{SymFourthOrderTensorValue{D,T,L}}) where {D,T,L} = length(SymFourthOrderTensorValue{D,T,L})
n_components(::SymFourthOrderTensorValue{D,T,L}) where {D,T,L} = n_components(SymFourthOrderTensorValue{D,T,L})

