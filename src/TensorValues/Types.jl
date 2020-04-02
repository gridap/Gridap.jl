###############################################################
# Types
###############################################################

"""
Type representing a multi-dimensional value
"""
abstract type MultiValue{S,T,N,L} <: Number end

"""
Type representing a first-order tensor
"""
struct VectorValue{D,T} <: MultiValue{Tuple{D},T,1,D}
    data::NTuple{D,T}
    function VectorValue{D,T}(data::NTuple{D,T}) where {D,T}
        new{D,T}(data)
    end
end

"""
Type representing a second-order tensor
"""
struct TensorValue{D1,D2,T,L} <: MultiValue{Tuple{D1,D2},T,2,L}
    data::NTuple{L,T}
    function TensorValue{D1,D2,T}(data::NTuple{L,T}) where {D1,D2,T,L}
        @assert L == D1*D2
        new{D1,D2,T,L}(data)
    end
end

###############################################################
# Constructors (VectorValue)
###############################################################

# Empty VectorValue constructor

VectorValue()                   = VectorValue{0,Int}(NTuple{0,Int}())
VectorValue{0}()                = VectorValue{0,Int}(NTuple{0,Int}())
VectorValue{0,T}() where {T}    = VectorValue{0,T}(NTuple{0,T}())
VectorValue(data::NTuple{0})    = VectorValue{0,Int}(data)
VectorValue{0}(data::NTuple{0}) = VectorValue{0,Int}(data)

# VectorValue single NTuple argument constructor

VectorValue(data::NTuple{D,T})        where {D,T}     = VectorValue{D,T}(data)
VectorValue{D}(data::NTuple{D,T})     where {D,T}     = VectorValue{D,T}(data)
VectorValue{D,T1}(data::NTuple{D,T2}) where {D,T1,T2} = VectorValue{D,T1}(NTuple{D,T1}(data))

# VectorValue Vararg constructor

VectorValue(data::Real...)                  = VectorValue(NTuple{length(data)}(data))
VectorValue{D}(data::Real...)   where {D}   = VectorValue{D}(NTuple{D}(data))
VectorValue{D,T}(data::Real...) where {D,T} = VectorValue{D,T}(NTuple{D,T}(data))

# VectorValue single SVector, MVector and AbstractVector argument constructor

function VectorValue(data::
                Union{
                    SVector{D,T2},
                    MVector{D,T2},
                    AbstractArray{T2}
                }) where {D,T1,T2}
    PD = (@isdefined D)  ? D  : length(data)
    VectorValue{PD,T2}(NTuple{PD,T2}(data))
end

function VectorValue{D}(data::
                Union{
                    SVector{D,T2},
                    MVector{D,T2},
                    AbstractArray{T2}
                }) where {D,T1,T2}
    VectorValue{D,T2}(NTuple{D,T2}(data))
end

function VectorValue{D,T1}(data::
                Union{
                    SVector{D,T2},
                    MVector{D,T2},
                    AbstractArray{T2}
                }) where {D,T1,T2}
    VectorValue{D,T1}(NTuple{D,T1}(data))
end

###############################################################
# Constructors (TensorValue)
###############################################################

# Empty TensorValue constructor

TensorValue()                     = TensorValue{0,0,Int}(NTuple{0,Int}())
TensorValue{0,0}()                = TensorValue{0,0,Int}(NTuple{0,Int}())
TensorValue{0,0,T}() where {T}    = TensorValue{0,0,T}(NTuple{0,T}())
TensorValue(data::NTuple{0})      = TensorValue{0,0,Int}(data)
TensorValue{0,0}(data::NTuple{0}) = TensorValue{0,0,Int}(data)

# TensorValue single NTuple argument constructor

TensorValue(data::NTuple{L,T}) where {L,T}                        = (D=Int(sqrt(L));TensorValue{D,D,T}(data))
TensorValue{D}(data::NTuple{L,T}) where {D,L,T}                   = TensorValue{D,D,T}(data)
TensorValue{D1,D2}(data::NTuple{L,T}) where {D1,D2,L,T}           = TensorValue{D1,D2,T}(data)
TensorValue{D1,D2,T1}(data::NTuple{L,T2}) where {D1,D2,L,T1,T2}   = TensorValue{D1,D2,T1}(NTuple{L,T1}(data))
TensorValue{D1,D2,T1,L}(data::NTuple{L,T2}) where {D1,D2,L,T1,T2} = TensorValue{D1,D2,T1}(NTuple{L,T1}(data))

# TensorValue Vararg constructor

TensorValue(data::Real...)                          = TensorValue(NTuple{length(data)}(data))
TensorValue{D}(data::Real...) where {D}             = TensorValue{D,D}(NTuple{D*D}(data))
TensorValue{D1,D2}(data::Real...) where {D1,D2}     = TensorValue{D1,D2}(NTuple{D1*D2}(data))
TensorValue{D1,D2,T}(data::Real...) where {D1,D2,T} = TensorValue{D1,D2,T}(NTuple{D1*D2,T}(data))

# VectorValue single SVector, MVector, SMatrix, MMatrix and AbstractMatrix argument constructor

function TensorValue(data::
                Union{
                    SMatrix{D1,D2,T2,L},
                    MMatrix{D1,D2,T2,L},
                    AbstractMatrix{T2}
                }) where {D1,D2,T1,T2,L}
    PD1 = (@isdefined D1) ? D1 : size(data)[1]
    PD2 = (@isdefined D2) ? D2 : size(data)[2]
    PL  = (@isdefined L)  ? L  : length(data)
    TensorValue{PD1,PD2,T2}(NTuple{PL,T2}(data))
end

function TensorValue{D1,D2}(data::
                Union{
                    SMatrix{D1,D2,T2,L},
                    MMatrix{D1,D2,T2,L},
                    AbstractMatrix{T2}
                }) where {D1,D2,T1,T2,L}
    PT  = (@isdefined T1) ? T1 : T2
    PL  = (@isdefined L)  ? L  : length(data)
    TensorValue{D1,D2,PT}(NTuple{PL,PT}(data))
end

function TensorValue{D1,D2,T1}(data::
                Union{
                    SMatrix{D1,D2,T2,L},
                    MMatrix{D1,D2,T2,L},
                    AbstractMatrix{T2}
                }) where {D1,D2,T1,T2,L}
    PL  = (@isdefined L)  ? L  : length(data)
    TensorValue{D1,D2,T1}(NTuple{PL,T1}(data))
end

function TensorValue{D1,D2,T1,L}(data::
                Union{
                    SMatrix{D1,D2,T2,L},
                    MMatrix{D1,D2,T2,L},
                    AbstractMatrix{T2}
                }) where {D1,D2,T1,T2,L}
    TensorValue{D1,D2,T1}(NTuple{L,T1}(data))
end

###############################################################
# Conversions (VectorValue)
###############################################################

function convert(::Type{<:Union{VectorValue,VectorValue{D,T1}}}, 
                arg::
                    Union{
                        NTuple{D,T2},
                        SVector{D,T2},
                        MVector{D,T2},
                        AbstractArray{T2}
                    }) where {D,T1,T2}
    PT = (@isdefined T1) ? T1 : T2
    PD = (@isdefined D)  ? D  : length(arg)
    VectorValue{PD,PT}(NTuple{PD,PT}(arg))
end

function convert(::Type{<:Union{NTuple,NTuple{D,T1}}}, arg::VectorValue{D,T2}) where {D,T1,T2}
    PT = (@isdefined T1) ? T1 : T2
    NTuple{D,PT}(arg.data)
end

function convert(::Type{<:Union{SVector,SVector{D,T1}}}, arg::VectorValue{D,T2}) where {D,T1,T2}
    PT = (@isdefined T1) ? T1 : T2
    SVector{D,PT}(arg.data)
end

function convert(::Type{<:Union{MVector,MVector{D,T1}}}, arg::VectorValue{D,T2}) where {D,T1,T2}
    PT = (@isdefined T1) ? T1 : T2
    MVector{D,PT}(arg.data)
end

function convert(::Type{<:Union{VectorValue,VectorValue{D,T1}}}, arg::VectorValue{D,T2}) where {D,T1,T2}
    PT = (@isdefined T1) ? T1 : T2
    PT == T2 ? arg : convert(VectorValue{D,PT}, arg.data)
end

###############################################################
# Conversions (TensorValue)
###############################################################

function convert(TT::Type{<:Union{TensorValue,TensorValue{D1,D2,T1},TensorValue{D1,D2,T1,L}}}, 
                arg::
                    Union{
                        NTuple{L,T2},
                        SMatrix{D1,D2,T2,L},
                        MMatrix{D1,D2,T2,L}
                    }) where {D1,D2,T1,T2,L}
    PT = (@isdefined T1) ? T1 : T2
    TensorValue{D1,D2,PT}(arg)
end

function convert(T::Type{<:Union{NTuple,NTuple{L,T1}}}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    PT = (@isdefined T1) ? T1 : T2
    NTuple{L,PT}(arg.data)
end

function convert(::Type{<:Union{SMatrix,SMatrix{D1,D2,T1,L}}}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    PT = (@isdefined T1) ? T1 : T2
    SMatrix{D1,D2,PT,L}(arg.data)
end

function convert(::Type{<:Union{MMatrix,MMatrix{D1,D2,T1,L}}}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    PT = (@isdefined T1) ? T1 : T2
    MMatrix{D1,D2,PT,L}(arg.data)
end

function convert(::Type{<:SMatrix{D1,D2,T1}}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    SMatrix{D1,D2,T1,L}(arg.data)
end

function convert(::Type{<:MMatrix{D1,D2,T1}}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    MMatrix{D1,D2,T1,L}(arg.data)
end

function convert(::Type{<:Union{TensorValue,TensorValue{D1,D2},TensorValue{D1,D2,T1},TensorValue{D1,D2,T1,L}}}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    PT = (@isdefined T1) ? T1 : T2
    PT == T2 ? arg : convert(TensorValue{D1,D2,T1,L}, arg.data)
end

###############################################################
# Other constructors and conversions (VectorValue)
###############################################################

zero(::Type{<:VectorValue{D,T}}) where {D,T} = VectorValue{D,T}(NTuple{D,T}(zeros(T,D)))
zero(::VectorValue{D,T}) where {D,T} = zero(VectorValue{D,T})

one(::Type{<:VectorValue{D,T}}) where {D,T} = VectorValue{D,T}(NTuple{D,T}(ones(T,D)))
one(::VectorValue{D,T}) where {D,T} = one(VectorValue{D,T})

mutable(::Type{VectorValue{D,T}}) where {D,T} = MVector{D,T}
mutable(::VectorValue{D,T}) where {D,T} = mutable(VectorValue{D,T})

change_eltype(::Type{VectorValue{D}},::Type{T}) where {D,T} = VectorValue{D,T}
change_eltype(::Type{VectorValue{D,T1}},::Type{T2}) where {D,T1,T2} = VectorValue{D,T2}
change_eltype(::VectorValue{D,T1},::Type{T2}) where {D,T1,T2} = change_eltype(VectorValue{D,T1},T2)

SVector(arg::VectorValue{D,T}) where {D,T} = SVector{D,T}(arg.data)
SArray(arg::VectorValue{D,T}) where {D,T} = SVector(arg)
get_array(arg::T where {T<:VectorValue}) = convert(SVector,arg)

###############################################################
# Other constructors and conversions (TensorValue)
###############################################################

zero(::Type{<:TensorValue{D1,D2,T}}) where {D1,D2,T} = (D=D1*D2;TensorValue{D1,D2,T}(NTuple{D,T}(zeros(T,D))))
zero(::TensorValue{D1,D2,T}) where {D1,D2,T} = zero(TensorValue{D1,D2,T})

@generated function one(::Type{<:TensorValue{D1,D2,T}}) where {D1,D2,T}
  str = join(["$i==$j ? one(T) : zero(T), " for i in 1:D1 for j in 1:D2])
  Meta.parse("TensorValue{D1,D2,T}(($str))")
end
one(::TensorValue{D1,D2,T}) where {D1,D2,T} = one(TensorValue{D1,D2,T})

mutable(::Type{<:TensorValue{D1,D2,T}}) where {D1,D2,T} = MMatrix{D1,D2,T}
mutable(::TensorValue{D1,D2,T}) where {D1,D2,T} = mutable(TensorValue{D1,D2,T})

change_eltype(::Type{TensorValue{D1,D2,T1,L}},::Type{T2}) where {D1,D2,T1,T2,L} = TensorValue{D1,D2,T2,L}
change_eltype(::TensorValue{D1,D2,T1,L},::Type{T2}) where {D1,D2,T1,T2,L} = change_eltype(TensorValue{D1,D2,T1,L},T2)

SMatrix(arg::TensorValue{D1,D2,T,L}) where {D1,D2,T,L} = SMatrix{D1,D2,T,L}(arg.data)
SArray(arg::TensorValue{D1,D2,T,L}) where {D1,D2,T,L} = StaticArrays.SMatrix(arg)
get_array(arg::T where {T<:TensorValue}) = convert(SMatrix,arg)

@generated function diagonal_tensor(v::VectorValue{D,T}) where {D,T}
  s = ["zero(T), " for i in 1:(D*D)]
  for i in 1:D
    d = D*(i-1)+i
    s[d] = "v.data[$i],"
  end
  str = join(s)
  Meta.parse("TensorValue(($str))")
end

###############################################################
# Introspection (VectorValue)
###############################################################

eltype(::Type{<:VectorValue{D,T}}) where {D,T} = T
eltype(arg::VectorValue{D,T}) where {D,T} = eltype(VectorValue{D,T})

size(::Type{VectorValue{D}}) where {D} = (D,)
size(::Type{VectorValue{D,T}}) where {D,T} = (D,)
size(::VectorValue{D,T}) where {D,T}  = size(VectorValue{D,T})

length(::Type{VectorValue{D}}) where {D} = D
length(::Type{VectorValue{D,T}}) where {D,T} = D
length(::VectorValue{D,T}) where {D,T} = length(VectorValue{D,T})

n_components(::Type{VectorValue{D}}) where {D} = length(VectorValue{D})
n_components(::Type{VectorValue{D,T}}) where {D,T} = length(VectorValue{D,T})
n_components(::VectorValue{D,T}) where {D,T} = n_components(VectorValue{D,T})

###############################################################
# Introspection (TensorValue)
###############################################################

eltype(::Type{<:TensorValue{D1,D2,T}}) where {D1,D2,T} = T
eltype(::TensorValue{D1,D2,T}) where {D1,D2,T} = eltype(TensorValue{D1,D2,T})

size(::Type{TensorValue{D}}) where {D} = (D,D)
size(::Type{TensorValue{D1,D2}}) where {D1,D2} = (D1,D2)
size(::Type{TensorValue{D1,D2,T}}) where {D1,D2,T} = (D1,D2)
size(::Type{TensorValue{D1,D2,T,L}}) where {D1,D2,T,L} = (D1,D2)
size(::TensorValue{D1,D2,T}) where {D1,D2,T} = size(TensorValue{D1,D2,T})

length(::Type{TensorValue{D}}) where {D} = length(TensorValue{D,D})
length(::Type{TensorValue{D1,D2}}) where {D1,D2} = D1*D1
length(::Type{TensorValue{D1,D2,T,L}}) where {D1,D2,T,L} = L
length(::TensorValue{D1,D2,T,L}) where {D1,D2,T,L} = length(TensorValue{D1,D2,T,L})

n_components(::Type{TensorValue{D}}) where {D} = length(TensorValue{D,D})
n_components(::Type{TensorValue{D1,D2}}) where {D1,D2} = length(TensorValue{D1,D2})
n_components(::Type{TensorValue{D1,D2,T,L}}) where {D1,D2,T,L} = length(TensorValue{D1,D2,T,L})
n_components(::TensorValue{D1,D2,T,L}) where {D1,D2,T,L} = n_components(TensorValue{D1,D2,T,L})


###############################################################
# Other constructors and conversions implemented for more general types
###############################################################

change_eltype(::Type{<:Number},::Type{T}) where {T} = T
change_eltype(::Number,::Type{T2}) where {T2} = change_eltype(Number,T2)

n_components(::Type{<:Number}) = 1
n_components(::Number) = n_components(Number)

###############################################################
# Misc
###############################################################

# Misc operations on the type itself

@pure _s(s::Size{T}) where T = T

# Custom type printing

function show(io::IO,v::MultiValue)
  print(io,v.data)
end

function show(io::IO,::MIME"text/plain",v:: MultiValue)
  print(io,typeof(v))
  print(io,v.data)
end

@inline Tuple(arg::MultiValue) = arg.data

