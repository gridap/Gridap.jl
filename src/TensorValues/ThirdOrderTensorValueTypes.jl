
"""
Type representing a third-order tensor
"""
struct ThirdOrderTensorValue{D1,D2,D3,T,L} <: MultiValue{Tuple{D1,D2,D3},T,3,L}
    data::NTuple{L,T}
    function ThirdOrderTensorValue{D1,D2,D3,T}(data::NTuple{L,T}) where {D1,D2,D3,T,L}
        @assert L == D1*D2*D3
        new{D1,D2,D3,T,L}(data)
    end
end

# Empty ThirdOrderTensorValue constructor

ThirdOrderTensorValue()                       = ThirdOrderTensorValue{0,0,0,Int}(NTuple{0,Int}())
ThirdOrderTensorValue{0,0,0}()                = ThirdOrderTensorValue{0,0,0,Int}(NTuple{0,Int}())
ThirdOrderTensorValue{0,0,0,T}() where{T}     = ThirdOrderTensorValue{0,0,0,T}(NTuple{0,T}())
ThirdOrderTensorValue(data::NTuple{0})        = ThirdOrderTensorValue{0,0,0,Int}(data)
ThirdOrderTensorValue{0,0,0}(data::NTuple{0}) = ThirdOrderTensorValue{0,0,0,Int}(data)

# ThirdOrderTensorValue single NTuple argument constructor

@generated function ThirdOrderTensorValue(data::NTuple{L,T}) where {L,T}
  D=Int(cbrt(L))
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

ThirdOrderTensorValue(data...) = ThirdOrderTensorValue(data)
ThirdOrderTensorValue{D}(data...) where {D} = ThirdOrderTensorValue{D}(data)
ThirdOrderTensorValue{D1,D2,D3}(data...) where {D1,D2,D3} = ThirdOrderTensorValue{D1,D2,D3}(data)
ThirdOrderTensorValue{D1,D2,D3,T1}(data...) where {D1,D2,D3,T1} =  ThirdOrderTensorValue{D1,D2,D3,T1}(data)
ThirdOrderTensorValue{D1,D2,D3,T1,L}(data...) where {D1,D2,D3,T1,L} =  ThirdOrderTensorValue{D1,D2,D3,T1}(data)

# From an array

ThirdOrderTensorValue(data::AbstractArray{T,3}) where {T} = ThirdOrderTensorValue(NTuple{length(data),T}(data))
ThirdOrderTensorValue{D}(data::AbstractArray{T,3}) where {D,T} = ThirdOrderTensorValue{D}(NTuple{length(data),T}(data))
ThirdOrderTensorValue{D,T1}(data::AbstractArray{T2,3}) where {D,T1,T2} = ThirdOrderTensorValue{D,T1}(NTuple{length(data),T}(data))

change_eltype(::Type{ThirdOrderTensorValue{D1,D2,D3,T1,L}},::Type{T2}) where {D1,D2,D3,T1,T2,L} = ThirdOrderTensorValue{D1,D2,D3,T2,L}
change_eltype(::T,::Type{T2}) where {T<:ThirdOrderTensorValue,T2} = change_eltype(T,T2)

zero(::Type{<:ThirdOrderTensorValue{D1,D2,D3,T}}) where {D1,D2,D3,T} = ThirdOrderTensorValue{D1,D2,D3,T}(tfill(zero(T),Val{D1*D2*D3}()))
zero(::ThirdOrderTensorValue{D1,D2,D3,T}) where {D1,D2,D3,T} = zero(ThirdOrderTensorValue{D1,D2,D3,T})
