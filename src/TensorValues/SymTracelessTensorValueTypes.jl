###############################################################
# SymTracelessTensorValue Type
###############################################################

"""
    SymTracelessTensorValue{D,T,L} <: AbstractSymTensorValue{D,T,L}
    QTensorValue{D,T,L}

Type representing a symetric second-order `D`×`D` tensor with zero trace. It must hold `L` = `D`(`D`+1)/2.
This type is used to model the Q-tensor order parameter in nematic liquid cristals.

The constructor determines the value of index (`D`,`D`) as minus the sum of the other diagonal values, so it value musn't be provided. The constructor thus expects the `L`-1 components of indices (i,j) for 1 ≤ i ≤ `D`-1 and  i ≤ j ≤ `D`.
"""
struct SymTracelessTensorValue{D,T,L} <: AbstractSymTensorValue{D,T,L}
  data::NTuple{L,T}
  function SymTracelessTensorValue{D,T}(data::NTuple{L,T}) where {D,T,L}
    @check L == D*(D+1)÷2-1
    new{D,T,L+1}( (data..., _minus_trace(data,Val(D))) )
  end
  function SymTracelessTensorValue{0,T}(data::NTuple{0,T}) where {T}
    new{0,T,0}(data)
  end
  function SymTracelessTensorValue{1,T}(::NTuple{0,T}) where {T}
    new{1,T,1}( (zero(T),) )
  end
end

@generated function _minus_trace(data::NTuple{L,T},::Val{D}) where {D,T,L}
  str = ""
  for i in 1:D-1
    k = _2d_sym_tensor_linear_index(D,i,i)
    str *= "- data[$k]"
  end
  Meta.parse("($str)")
end

const QTensorValue = SymTracelessTensorValue

function promote_rule(
  ::Type{<:SymTracelessTensorValue{D,Ta}},
  ::Type{<:SymTracelessTensorValue{D,Tb}}) where {D,Ta,Tb}

  T = promote_type(Ta,Tb)
  SymTracelessTensorValue{D,T}
end

###############################################################
# Constructors (SymTracelessTensorValue)
###############################################################

# Empty SymTracelessTensorValue constructor

SymTracelessTensorValue()                   = SymTracelessTensorValue{0,Int}(NTuple{0,Int}())
SymTracelessTensorValue{0}()                = SymTracelessTensorValue{0,Int}(NTuple{0,Int}())
SymTracelessTensorValue{0,T}() where {T}    = SymTracelessTensorValue{0,T}(NTuple{0,T}())
SymTracelessTensorValue(data::NTuple{0})    = SymTracelessTensorValue{0,Int}(data)
SymTracelessTensorValue{0}(data::NTuple{0}) = SymTracelessTensorValue{0,Int}(data)

# 1D SymTracelessTensorValue missing constructor

SymTracelessTensorValue{1}()                = SymTracelessTensorValue{1,Int}(NTuple{0,Int}())
SymTracelessTensorValue{1}(data::NTuple{0}) = SymTracelessTensorValue{1,Int}(data)

# SymTracelessTensorValue single NTuple argument constructor

@generated function SymTracelessTensorValue(data::NTuple{L,T}) where {L,T}
  msg = "Invalid number of scalar arguments in SymTracelessTensorValue constructor"
  V = (sqrt(9+8*L)-1)/2
  @check floor(Int,V) == ceil(Int,V) msg
  D = Int(V)
  quote
    SymTracelessTensorValue{$D,T}(data)
  end
end
SymTracelessTensorValue{D}(data::NTuple{L,T}) where {D,L,T}           = SymTracelessTensorValue{D,T}(data)
SymTracelessTensorValue{D,T1}(data::NTuple{L,T2}) where {D,L,T1,T2}   = SymTracelessTensorValue{D,T1}(NTuple{L,T1}(data))
SymTracelessTensorValue{D,T1,L}(data::NTuple{Lm,T2}) where {D,L,Lm,T1,T2} = SymTracelessTensorValue{D,T1}(NTuple{Lm,T1}(data))

# SymTracelessTensorValue single Tuple argument constructor

SymTracelessTensorValue(data::Tuple)                      = SymTracelessTensorValue(promote(data...))
SymTracelessTensorValue{D}(data::Tuple) where {D}         = SymTracelessTensorValue{D}(promote(data...))
SymTracelessTensorValue{D,T1}(data::Tuple) where {D,T1}   = SymTracelessTensorValue{D,T1}(NTuple{length(data),T1}(data))
SymTracelessTensorValue{D,T1,L}(data::Tuple) where {D,T1,L} = SymTracelessTensorValue{D,T1}(NTuple{L-1,T1}(data))

# SymTracelessTensorValue Vararg constructor

SymTracelessTensorValue(data::Number...)                    = SymTracelessTensorValue(data)
SymTracelessTensorValue{D}(data::Number...)    where {D}    = SymTracelessTensorValue{D}(data)
SymTracelessTensorValue{D,T1}(data::Number...) where {D,T1} = SymTracelessTensorValue{D,T1}(data)
SymTracelessTensorValue{D,T1,L}(data::Number...) where {D,T1,L} = SymTracelessTensorValue{D,T1}(data)

# SymTracelessTensorValue single AbstractMatrix argument constructor

#From Square Matrices
@generated function _flatten_upper_triangle_traceless(data::AbstractArray,::Val{D}) where D
  check_e = :( @check ($D,$D) == size(data) )
  str = ""
  for i in 1:D-1
    for j in i:D
      str *= "data[$i,$j], "
    end
  end
  ret_e = Meta.parse(" return ($str)")
  Expr(:block, check_e, ret_e)
end

SymTracelessTensorValue(data::AbstractMatrix{T}) where {T} = ((D1,D2)=size(data); SymTracelessTensorValue{D1}(data))
SymTracelessTensorValue{D}(data::AbstractMatrix{T}) where {D,T} = SymTracelessTensorValue{D,T}(_flatten_upper_triangle_traceless(data,Val{D}()))
SymTracelessTensorValue{D,T1}(data::AbstractMatrix{T2}) where {D,T1,T2} = SymTracelessTensorValue{D,T1}(_flatten_upper_triangle_traceless(data,Val{D}()))
SymTracelessTensorValue{D,T1,L}(data::AbstractMatrix{T2}) where {D,T1,T2,L} = SymTracelessTensorValue{D,T1,L}(_flatten_upper_triangle_traceless(data,Val{D}()))

###############################################################
# Conversions (SymTracelessTensorValue)
###############################################################

@generated function _SymTracelessTensorValue_to_array(arg::SymTracelessTensorValue{D,T,L}) where {D,T,L}
  str = ""
  for j in 1:D
    for i in 1:D
      str *= "arg[$i,$j], "
    end
  end
  Meta.parse("SMatrix{D,D,T}(($str))")
end

# Inverse conversion
convert(::Type{<:MArray{Tuple{D,D},T}}, arg::SymTracelessTensorValue) where {D,T} = MMatrix{D,D,T}(_SymTracelessTensorValue_to_array(arg))
convert(::Type{<:SArray{Tuple{D,D},T}}, arg::SymTracelessTensorValue) where {D,T} = _SymTracelessTensorValue_to_array(arg)

# Internal conversion
convert(::Type{<:SymTracelessTensorValue{D}}, arg::SymTracelessTensorValue{D}) where {D} = SymTracelessTensorValue{D}(Tuple(arg)[1:end-1])
convert(::Type{<:SymTracelessTensorValue{D,T}}, arg::SymTracelessTensorValue{D}) where {D,T} = SymTracelessTensorValue{D,T}(Tuple(arg)[1:end-1])

###############################################################
# Other constructors and conversions (SymTracelessTensorValue)
###############################################################

change_eltype(::Type{<:SymTracelessTensorValue{D}},::Type{T2}) where {D,T2} = SymTracelessTensorValue{D,T2}
change_eltype(::Type{SymTracelessTensorValue{D,T1,L}},::Type{T2}) where {D,T1,T2,L} = SymTracelessTensorValue{D,T2,L}

num_indep_components(::Type{<:SymTracelessTensorValue{0}}) = 0
num_indep_components(::Type{<:SymTracelessTensorValue{D}}) where {D} = D*(D+1)÷2-1

