###############################################################
# Types
###############################################################

"""
Type representing a first-order tensor
"""
struct VectorValue{D,T} <: MultiValue{Tuple{D},T,1,D}
    data::NTuple{D,T}
    function VectorValue{D,T}(data::NTuple{D,T}) where {D,T}
        new{D,T}(data)
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

VectorValue{D}(data::NTuple{D2,T2}) where {D,D2,T2} = @unreachable
VectorValue{D1,T1}(data::NTuple{D2,T2}) where {D1,T1,D2,T2} = @unreachable

# VectorValue single Tuple argument constructor

VectorValue(data::Tuple)                    = VectorValue(promote(data...))
VectorValue{D}(data::Tuple)    where {D}    = VectorValue{D}(promote(data...))
VectorValue{D,T1}(data::Tuple) where {D,T1} = VectorValue{D,T1}(NTuple{D,T1}(data))

# VectorValue Vararg constructor

VectorValue(data...)                    = VectorValue(data)
VectorValue{D}(data...)    where {D}    = VectorValue{D}(data)
VectorValue{D,T1}(data...) where {D,T1} = VectorValue{D,T1}(data)

# Fix for julia 1.0.4
VectorValue{D}(data::T...)    where {D,T}    = VectorValue{D,T}(data)

# VectorValue single AbstractVector argument constructor

VectorValue(data::AbstractArray{T}) where {T}              = (D=length(data);VectorValue(NTuple{D,T}(data)))
VectorValue{D}(data::AbstractArray{T}) where {D,T}         = VectorValue{D}(NTuple{D,T}(data))
VectorValue{D,T1}(data::AbstractArray{T2}) where {D,T1,T2} = VectorValue{D,T1}(NTuple{D,T1}(data))

###############################################################
# Conversions (VectorValue)
###############################################################

# Direct conversion
convert(::Type{<:VectorValue{D,T}}, arg:: AbstractArray) where {D,T} = VectorValue{D,T}(NTuple{D,T}(arg))
convert(::Type{<:VectorValue{D,T}}, arg:: Tuple) where {D,T} = VectorValue{D,T}(arg)

# Inverse conversion
convert(::Type{<:SVector{D,T}}, arg::VectorValue{D}) where {D,T} = SVector{D,T}(Tuple(arg))
convert(::Type{<:MVector{D,T}}, arg::VectorValue{D}) where {D,T} = MVector{D,T}(Tuple(arg))
convert(::Type{<:NTuple{D,T}},  arg::VectorValue{D}) where {D,T} = NTuple{D,T}(Tuple(arg))

# Internal conversion
convert(::Type{<:VectorValue{D,T}}, arg::VectorValue{D}) where {D,T} = VectorValue{D,T}(Tuple(arg))
convert(::Type{<:VectorValue{D,T}}, arg::VectorValue{D,T}) where {D,T} = arg

###############################################################
# Other constructors and conversions (VectorValue)
###############################################################

zero(::Type{<:VectorValue{D,T}}) where {D,T} = VectorValue{D,T}(tfill(zero(T),Val{D}()))
zero(::VectorValue{D,T}) where {D,T} = zero(VectorValue{D,T})

Mutable(::Type{VectorValue{D,T}}) where {D,T} = MVector{D,T}
Mutable(::VectorValue{D,T}) where {D,T} = Mutable(VectorValue{D,T})
mutable(a::VectorValue) = MVector(a.data)

change_eltype(::Type{VectorValue{D}},::Type{T}) where {D,T} = VectorValue{D,T}
change_eltype(::Type{VectorValue{D,T1}},::Type{T2}) where {D,T1,T2} = VectorValue{D,T2}
change_eltype(::VectorValue{D,T1},::Type{T2}) where {D,T1,T2} = change_eltype(VectorValue{D,T1},T2)

get_array(arg::VectorValue{D,T}) where {D,T} = convert(SVector{D,T}, arg)

###############################################################
# Introspection (VectorValue)
###############################################################

eltype(::Type{<:VectorValue{D,T}}) where {D,T} = T
eltype(arg::VectorValue{D,T}) where {D,T} = eltype(VectorValue{D,T})

size(::Type{<:VectorValue{D}}) where {D} = (D,)
size(::VectorValue{D}) where {D}  = size(VectorValue{D})

length(::Type{<:VectorValue{D}}) where {D} = D
length(::VectorValue{D}) where {D} = length(VectorValue{D})

num_components(::Type{<:VectorValue{D}}) where {D} = length(VectorValue{D})
num_components(::VectorValue{D}) where {D} = num_components(VectorValue{D})

