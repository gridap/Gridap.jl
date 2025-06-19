###############################################################
# SkewSymTensorValue Type
###############################################################

"""
    SkewSymTensorValue{D,T,L} <: MultiValue{Tuple{D,D},T,2,L}

Type representing a skew symmetric second-order `D`×`D` tensor.
It must hold `L` = `D`(`D`-1)/2.

It is constructed by providing the components of index (i,j) for 1 ≤ i < j ≤ `D`.
"""
struct SkewSymTensorValue{D,T,L} <: MultiValue{Tuple{D,D},T,2,L}
  data::NTuple{L,T}
  function SkewSymTensorValue{D,T}(data::NTuple{L,T}) where {D,T,L}
    @check L == D*(D-1)÷2
    new{D,T,L}(data)
  end
end

function promote_rule(::Type{<:SkewSymTensorValue{D,Ta}}, ::Type{<:SkewSymTensorValue{D,Tb}}) where {D,Ta,Tb}
  T = promote_type(Ta,Tb)
  SkewSymTensorValue{D,T}
end

###############################################################
# Constructors (SkewSymTensorValue)
###############################################################

# Empty SkewSymTensorValue constructor

SkewSymTensorValue()                   = SkewSymTensorValue{0,Int}(NTuple{0,Int}())
SkewSymTensorValue{0}()                = SkewSymTensorValue{0,Int}(NTuple{0,Int}())
SkewSymTensorValue{0,T}() where {T}    = SkewSymTensorValue{0,T}(NTuple{0,T}())
SkewSymTensorValue(data::NTuple{0})    = SkewSymTensorValue{0,Int}(data)
SkewSymTensorValue{0}(data::NTuple{0}) = SkewSymTensorValue{0,Int}(data)

# SkewSymTensorValue single NTuple argument constructor

@generated function SkewSymTensorValue(data::NTuple{L,T}) where {L,T}
  msg = "Invalid number of scalar arguments in SkewSymTensorValue constructor"
  V = (sqrt(1+8*L)+1)/2
  @check floor(Int,V) == ceil(Int,V) msg
  D = Int(V)
  quote
    SkewSymTensorValue{$D,T}(data)
  end
end
SkewSymTensorValue{D}(data::NTuple{L,T}) where {D,L,T}           = SkewSymTensorValue{D,T}(data)
SkewSymTensorValue{D,T1}(data::NTuple{L,T2}) where {D,L,T1,T2}   = SkewSymTensorValue{D,T1}(NTuple{L,T1}(data))
SkewSymTensorValue{D,T1,L}(data::NTuple{L,T2}) where {D,L,T1,T2} = SkewSymTensorValue{D,T1}(NTuple{L,T1}(data))

# SkewSymTensorValue single Tuple argument constructor

SkewSymTensorValue(data::Tuple)                      = SkewSymTensorValue(promote(data...))
SkewSymTensorValue{D}(data::Tuple) where {D}         = SkewSymTensorValue{D}(promote(data...))
SkewSymTensorValue{D,T1}(data::Tuple) where {D,T1}   = SkewSymTensorValue{D,T1}(NTuple{length(data),T1}(data))
SkewSymTensorValue{D,T1,L}(data::Tuple) where {D,T1,L}   = SkewSymTensorValue{D,T1}(NTuple{L,T1}(data))

# SkewSymTensorValue Vararg constructor

SkewSymTensorValue(data::Number...)                    = SkewSymTensorValue(data)
SkewSymTensorValue{D}(data::Number...)    where {D}    = SkewSymTensorValue{D}(data)
SkewSymTensorValue{D,T1}(data::Number...) where {D,T1} = SkewSymTensorValue{D,T1}(data)
SkewSymTensorValue{D,T1,L}(data::Number...) where {D,T1,L} = SkewSymTensorValue{D,T1}(data)

# SkewSymTensorValue single AbstractMatrix argument constructor

#From Square Matrices
@generated function _flatten_skewsym(data::AbstractArray,::Val{D}) where D
  check_e = :( @check ($D,$D) == size(data) )
  str = ""
  for i in 1:D
    for j in i+1:D
      str *= "data[$i,$j], "
    end
  end
  ret_e = Meta.parse(" return ($str)")
  Expr(:block, check_e, ret_e)
end

SkewSymTensorValue(data::AbstractMatrix{T}) where {T} = ((D1,D2)=size(data); SkewSymTensorValue{D1}(data))
SkewSymTensorValue{D}(data::AbstractMatrix{T}) where {D,T} = SkewSymTensorValue{D,T}(_flatten_skewsym(data,Val{D}()))
SkewSymTensorValue{D,T1}(data::AbstractMatrix{T2}) where {D,T1,T2} = SkewSymTensorValue{D,T1}(_flatten_skewsym(data,Val{D}()))
SkewSymTensorValue{D,T1,L}(data::AbstractMatrix{T2}) where {D,T1,T2,L} = SkewSymTensorValue{D,T1,L}(_flatten_skewsym(data,Val{D}()))

###############################################################
# Conversions (SkewSymTensorValue)
###############################################################

@generated function _SkewSymTensorValue_to_array(arg::SkewSymTensorValue{D,T,L}) where {D,T,L}
  z = zero(T)
  str = ""
  for j in 1:D
    for i in 1:D
      str *= "arg[$i,$j], "
    end
  end
  Meta.parse("SMatrix{D,D,T}(($str))")
end

# Inverse conversion
convert(::Type{<:MArray{Tuple{D,D},T}}, arg::SkewSymTensorValue) where {D,T} = MMatrix{D,D,T}(_SkewSymTensorValue_to_array(arg))
convert(::Type{<:SArray{Tuple{D,D},T}}, arg::SkewSymTensorValue) where {D,T} = _SkewSymTensorValue_to_array(arg)

# Internal conversion
convert(::Type{<:SkewSymTensorValue{D,T}}, arg::SkewSymTensorValue{D}) where {D,T} = SkewSymTensorValue{D,T}(Tuple(arg))
convert(::Type{<:SkewSymTensorValue{D,T}}, arg::SkewSymTensorValue{D,T}) where {D,T} = arg

###############################################################
# Other constructors and conversions (SkewSymTensorValue)
###############################################################

one(::Type{<:SkewSymTensorValue}) = @unreachable "Skew symmetric tensors do not have multiplicative neutral."

change_eltype(::Type{<:SkewSymTensorValue{D}},::Type{T2}) where {D,T2} = SkewSymTensorValue{D,T2}
change_eltype(::Type{SkewSymTensorValue{D,T1,L}},::Type{T2}) where {D,T1,T2,L} = SkewSymTensorValue{D,T2,L}

num_indep_components(::Type{<:SkewSymTensorValue{D}}) where {D} = D*(D-1)÷2

