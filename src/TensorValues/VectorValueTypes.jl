###############################################################
# Types
###############################################################

"""
    VectorValue{D,T} <: MultiValue{Tuple{D},T,1,D}

Type representing a first-order tensor, that is a vector, of length `D`.
"""
struct VectorValue{D,T} <: MultiValue{Tuple{D},T,1,D}
    data::NTuple{D,T}
    function VectorValue{D,T}(data::NTuple{D,T}) where {D,T}
        new{D,T}(data)
    end
end

function promote_rule(::Type{VectorValue{D,Ta}}, ::Type{VectorValue{D,Tb}}) where {D,Ta,Tb}
  VectorValue{D,promote_type(Ta,Tb)}
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

VectorValue(data::Number...)                    = VectorValue(data)
VectorValue{D}(data::Number...)    where {D}    = VectorValue{D}(data)
VectorValue{D,T1}(data::Number...) where {D,T1} = VectorValue{D,T1}(data)

# Fix for julia 1.0.4
VectorValue{D}(data::T...)    where {D,T<:Number}    = VectorValue{D,T}(data)

# VectorValue single AbstractVector argument constructor

VectorValue(data::AbstractArray{T}) where {T}              = (D=length(data);VectorValue(NTuple{D,T}(data)))
VectorValue{D}(data::AbstractArray{T}) where {D,T}         = VectorValue{D}(NTuple{D,T}(data))
VectorValue{D,T1}(data::AbstractArray{T2}) where {D,T1,T2} = VectorValue{D,T1}(NTuple{D,T1}(data))

###############################################################
# Conversions (VectorValue)
###############################################################

# Inverse conversion
convert(::Type{<:SArray{Tuple{D},T}}, arg::VectorValue{D}) where {D,T} = SVector{D,T}(Tuple(arg))
convert(::Type{<:MArray{Tuple{D},T}}, arg::VectorValue{D}) where {D,T} = MVector{D,T}(Tuple(arg))

# Internal conversion
convert(::Type{<:VectorValue{D,T}}, arg::VectorValue{D}) where {D,T} = VectorValue{D,T}(Tuple(arg))
convert(::Type{<:VectorValue{D,T}}, arg::VectorValue{D,T}) where {D,T} = arg

change_eltype(::Type{<:VectorValue{D}},::Type{T}) where {D,T} = VectorValue{D,T}

# Construction from SArray or MArray
MultiValue(a::StaticVector{D,T}) where {D,T} = convert(VectorValue{D,T}, a)

###############################################################
# VTK export (VectorValue)
###############################################################

function indep_components_names(::Type{<:VectorValue{A}}) where A
  [ "$i" for i in 1:A ]
  if A>3
    return ["$i" for i in 1:A ]
  else
    c_name = ["X", "Y", "Z"]
    return [c_name[i] for i in 1:A ]
  end
end
