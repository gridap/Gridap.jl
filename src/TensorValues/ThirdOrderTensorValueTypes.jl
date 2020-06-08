
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

# ThirdOrderTensorValue Vararg constructor

ThirdOrderTensorValue(data::T...) where {T}                      = ThirdOrderTensorValue(data)
ThirdOrderTensorValue{D}(data::T...) where {D,T}                 = ThirdOrderTensorValue{D}(data)
ThirdOrderTensorValue{D1,D2,D3}(data::T...) where {D1,D2,D3,T}   = ThirdOrderTensorValue{D1,D2,D3}(data)
ThirdOrderTensorValue{D1,D2,D3,T1}(data::T2...) where {D1,D2,D3,T1,T2} =  ThirdOrderTensorValue{D1,D2,D3,T1}(data)
ThirdOrderTensorValue{D1,D2,D3,T1,L}(data::T2...) where {D1,D2,D3,L,T1,T2} = ThirdOrderTensorValue{D1,D2,D3,T1,L}(data)

change_eltype(::Type{ThirdOrderTensorValue{D1,D2,D3,T1,L}},::Type{T2}) where {D1,D2,D3,T1,T2,L} = ThirdOrderTensorValue{D1,D2,D3,T2,L}
change_eltype(::T,::Type{T2}) where T<:ThirdOrderTensorValue = change_eltype(T,T2)
