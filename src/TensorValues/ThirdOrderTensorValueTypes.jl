
"""
    ThirdOrderTensorValue{D1,D2,D3,T,L} <: MultiValue{Tuple{D1,D2,D3},T,3,L}

Type representing a third-order `D1`×`D2`×`D3` tensor. It must hold `L` = `D1`\\*`D2`\\*`D3`.

If only `D1` or no dimension parameter is given to the constructor, `D1`=`D2`=`D3` is assumed.
"""
struct ThirdOrderTensorValue{D1,D2,D3,T,L} <: MultiValue{Tuple{D1,D2,D3},T,3,L}
    data::NTuple{L,T}
    function ThirdOrderTensorValue{D1,D2,D3,T}(data::NTuple{L,T}) where {D1,D2,D3,T,L}
        @check L == D1*D2*D3
        new{D1,D2,D3,T,L}(data)
    end
end

function promote_rule(
  ::Type{<:ThirdOrderTensorValue{D1,D2,D3,Ta}},
  ::Type{<:ThirdOrderTensorValue{D1,D2,D3,Tb}}) where {D1,D2,D3,Ta,Tb}

  T = promote_type(Ta,Tb)
  ThirdOrderTensorValue{D1,D2,D3,T}
end

# Empty ThirdOrderTensorValue constructor

ThirdOrderTensorValue()                       = ThirdOrderTensorValue{0,0,0,Int}(NTuple{0,Int}())
ThirdOrderTensorValue{0,0,0}()                = ThirdOrderTensorValue{0,0,0,Int}(NTuple{0,Int}())
ThirdOrderTensorValue{0,0,0,T}() where{T}     = ThirdOrderTensorValue{0,0,0,T}(NTuple{0,T}())
ThirdOrderTensorValue(data::NTuple{0})        = ThirdOrderTensorValue{0,0,0,Int}(data)
ThirdOrderTensorValue{0,0,0}(data::NTuple{0}) = ThirdOrderTensorValue{0,0,0,Int}(data)

# ThirdOrderTensorValue single NTuple argument constructor

@generated function ThirdOrderTensorValue(data::NTuple{L,T}) where {L,T}
  msg = "Invalid number of scalar arguments in ThirdOrderTensorValue constructor, expecting a cube number, got L=$L"
  V = cbrt(L)
  @assert floor(Int,V) == ceil(Int,V) msg
  D=Int(V)
  quote
    ThirdOrderTensorValue{$D,$D,$D,T}(data)
  end
end
ThirdOrderTensorValue{D}(data::NTuple{L,T}) where {D,L,T}                         = ThirdOrderTensorValue{D,D,D,T}(data)
ThirdOrderTensorValue{D1,D2,D3}(data::NTuple{L,T}) where {D1,D2,D3,L,T}           = ThirdOrderTensorValue{D1,D2,D3,T}(data)
ThirdOrderTensorValue{D1,D2,D3,T1}(data::NTuple{L,T2}) where {D1,D2,D3,L,T1,T2}   = ThirdOrderTensorValue{D1,D2,D3,T1}(NTuple{L,T1}(data))
ThirdOrderTensorValue{D1,D2,D3,T1,L}(data::NTuple{L,T2}) where {D1,D2,D3,L,T1,T2} = ThirdOrderTensorValue{D1,D2,D3,T1}(NTuple{L,T1}(data))

# ThirdOrderTensorValue single Tuple argument constructor

ThirdOrderTensorValue(data::Tuple) = ThirdOrderTensorValue(promote(data...))
ThirdOrderTensorValue{D}(data::Tuple) where {D} = ThirdOrderTensorValue{D,D,D}(promote(data...))
ThirdOrderTensorValue{D1,D2,D3}(data::Tuple) where {D1,D2,D3} = ThirdOrderTensorValue{D1,D2,D3}(promote(data...))
ThirdOrderTensorValue{D1,D2,D3,T1}(data::Tuple) where {D1,D2,D3,T1}  = ThirdOrderTensorValue{D1,D2,D3,T1}(NTuple{length(data),T1}(data))
ThirdOrderTensorValue{D1,D2,D3,T1,L}(data::Tuple) where {D1,D2,D3,T1,L}  = ThirdOrderTensorValue{D1,D2,D3,T1}(NTuple{L,T1}(data))

# ThirdOrderTensorValue Vararg constructor

ThirdOrderTensorValue(data::Number...) = ThirdOrderTensorValue(data)
ThirdOrderTensorValue{D}(data::Number...) where {D} = ThirdOrderTensorValue{D}(data)
ThirdOrderTensorValue{D1,D2,D3}(data::Number...) where {D1,D2,D3} = ThirdOrderTensorValue{D1,D2,D3}(data)
ThirdOrderTensorValue{D1,D2,D3,T1}(data::Number...) where {D1,D2,D3,T1} =  ThirdOrderTensorValue{D1,D2,D3,T1}(data)
ThirdOrderTensorValue{D1,D2,D3,T1,L}(data::Number...) where {D1,D2,D3,T1,L} =  ThirdOrderTensorValue{D1,D2,D3,T1}(data)

# ThirdOrderTensorValue single AbstractArray{3,T} argument constructor

ThirdOrderTensorValue(data::AbstractArray{T,3}) where {T} = ((D1,D2,D3)=size(data);L=length(data);ThirdOrderTensorValue{D1,D2,D3,T}(NTuple{L,T}(data)))
ThirdOrderTensorValue{D}(data::AbstractArray{T,3}) where {D,T} = (L=length(data);ThirdOrderTensorValue{D,D,D,T}(NTuple{L,T}(data)))
ThirdOrderTensorValue{D1,D2,D3}(data::AbstractArray{T,3}) where {D1,D2,D3,T} = (L=length(data);ThirdOrderTensorValue{D1,D2,D3,T}(NTuple{L,T}(data)))
ThirdOrderTensorValue{D1,D2,D3,T1}(data::AbstractArray{T2,3}) where {D1,D2,D3,T1,T2} = (L=length(data);ThirdOrderTensorValue{D1,D2,D3,T1}(NTuple{L,T1}(data)))
ThirdOrderTensorValue{D1,D2,D3,T1,L}(data::AbstractArray{T2,3}) where {D1,D2,D3,T1,T2,L} = ThirdOrderTensorValue{D1,D2,D3,T1}(NTuple{L,T1}(data))

###############################################################
# Conversions (ThirdOrderTensorValue)
###############################################################

# Inverse conversion
convert(::Type{<:SArray{Tuple{D1,D2,D3},T}}, arg::ThirdOrderTensorValue) where {D1,D2,D3,T} = SArray{Tuple{D1,D2,D3},T}(Tuple(arg))
convert(::Type{<:MArray{Tuple{D1,D2,D3},T}}, arg::ThirdOrderTensorValue) where {D1,D2,D3,T} = MArray{Tuple{D1,D2,D3},T}(Tuple(arg))

# Internal conversion
convert(::Type{<:ThirdOrderTensorValue{D1,D2,D3,T}}, arg::ThirdOrderTensorValue{D1,D2,D3}) where {D1,D2,D3,T} = ThirdOrderTensorValue{D1,D2,D3,T}(Tuple(arg))
convert(::Type{<:ThirdOrderTensorValue{D1,D2,D3,T}}, arg::ThirdOrderTensorValue{D1,D2,D3,T}) where {D1,D2,D3,T} = arg

# Construction from SArray or MArray
MultiValue(a::StaticArray{Tuple{D1,D2,D3},T}) where {D1,D2,D3,T} = convert(ThirdOrderTensorValue{D1,D2,D3,T}, a)

# other
change_eltype(::Type{<:ThirdOrderTensorValue{D1,D2,D3}},::Type{T2}) where {D1,D2,D3,T2} = ThirdOrderTensorValue{D1,D2,D3,T2}
change_eltype(::Type{ThirdOrderTensorValue{D1,D2,D3,T1,L}},::Type{T2}) where {D1,D2,D3,T1,T2,L} = ThirdOrderTensorValue{D1,D2,D3,T2,L}

###############################################################
# VTK export (ThirdOrderTensorValue)
###############################################################

function indep_components_names(::Type{<:ThirdOrderTensorValue{D1,D2,D3}}) where {D1,D2,D3}
  if D1>3 || D2>3 || D3>3
    return ["$i$j$k" for i in 1:D1 for j in 1:D2 for k in 1:D3 ]
  else
    c_name = ["X", "Y", "Z"]
    return [c_name[i]*c_name[j]*c_name[k] for i in 1:D1 for j in 1:D2 for k in 1:D3]
  end
end

