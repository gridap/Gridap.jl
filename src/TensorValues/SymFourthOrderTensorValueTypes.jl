###############################################################
# SymFourthOrderTensorValue Type
###############################################################

"""
    SymFourthOrderTensorValue{D,T,L} <: MultiValue{Tuple{D,D,D,D},T,4,L}

Type representing a symmetric second-order `D`×`D`×`D`×`D` tensor, with symmetries ijkl↔jikl and ijkl↔ijlk. It must hold `L` = (`D`(`D`+1)/2)^2.

It is constructed by providing the components of index (i,j,k,l) for 1 ≤ i ≤ j ≤ `D` and 1 ≤ k ≤ l ≤ `D`.
"""
struct SymFourthOrderTensorValue{D,T,L} <: MultiValue{Tuple{D,D,D,D},T,4,L}
  data::NTuple{L,T}
  function SymFourthOrderTensorValue{D,T}(data::NTuple{L,T}) where {D,T,L}
    @check L == (D*(D+1)÷2)^2
    new{D,T,L}(data)
  end
end

function promote_rule(::Type{<:SymFourthOrderTensorValue{D,Ta}}, ::Type{<:SymFourthOrderTensorValue{D,Tb}}) where {D,Ta,Tb}
  T = promote_type(Ta,Tb)
  SymFourthOrderTensorValue{D,T}
end

###############################################################
# Constructors (SymFourthOrderTensorValue)
###############################################################

# Empty SymFourthOrderTensorValue constructor

SymFourthOrderTensorValue()                   = SymFourthOrderTensorValue{0,Int}(NTuple{0,Int}())
SymFourthOrderTensorValue{0}()                = SymFourthOrderTensorValue{0,Int}(NTuple{0,Int}())
SymFourthOrderTensorValue{0,T}() where {T}    = SymFourthOrderTensorValue{0,T}(NTuple{0,T}())
SymFourthOrderTensorValue(data::NTuple{0})    = SymFourthOrderTensorValue{0,Int}(data)
SymFourthOrderTensorValue{0}(data::NTuple{0}) = SymFourthOrderTensorValue{0,Int}(data)

# SymFourthOrderTensorValue single NTuple argument constructor

@generated function SymFourthOrderTensorValue(data::NTuple{L,T}) where {L,T}
  msg = "Invalid number of scalar arguments in SymFourthOrderTensorValue constructor"
  V = (sqrt(1+8*sqrt(L))-1)/2
  @assert floor(Int,V) == ceil(Int,V) msg
  D = Int(V)
  quote
    SymFourthOrderTensorValue{$D,T}(data)
  end
end
SymFourthOrderTensorValue{D}(data::NTuple{L,T}) where {D,L,T}           = SymFourthOrderTensorValue{D,T}(data)
SymFourthOrderTensorValue{D,T1}(data::NTuple{L,T2}) where {D,L,T1,T2}   = SymFourthOrderTensorValue{D,T1}(NTuple{L,T1}(data))
SymFourthOrderTensorValue{D,T1,L}(data::NTuple{L,T2}) where {D,L,T1,T2} = SymFourthOrderTensorValue{D,T1}(NTuple{L,T1}(data))

# SymFourthOrderTensorValue single Tuple argument constructor

SymFourthOrderTensorValue(data::Tuple) = SymFourthOrderTensorValue(promote(data...))
SymFourthOrderTensorValue{D}(data::Tuple) where {D} = SymFourthOrderTensorValue{D}(promote(data...))
SymFourthOrderTensorValue{D,T1}(data::Tuple) where {D,T1} = SymFourthOrderTensorValue{D,T1}(NTuple{length(data),T1}(data))

# SymFourthOrderTensorValue Vararg constructor

SymFourthOrderTensorValue(data::Number...) = SymFourthOrderTensorValue(data)
SymFourthOrderTensorValue{D}(data::Number...) where {D} = SymFourthOrderTensorValue{D}(data)
SymFourthOrderTensorValue{D,T1}(data::Number...) where {D,T1} = SymFourthOrderTensorValue{D,T1}(data)

# SymFourthOrderTensorValue single AbstractArray argument constructor

#From square 4-dim Array
@generated function _flatten_sym_fourth_order_tensor(data::AbstractArray,::Val{D}) where D
  check_e = :( @check ($D,$D,$D,$D) == size(data) )
  str = ""
  for i in 1:D
    for j in i:D
      for k in 1:D
        for l in k:D
          str *= "data[$i,$j,$k,$l],"
        end
      end
    end
  end
  ret_e = Meta.parse("return ($str)")
  Expr(:block, check_e, ret_e)
end

SymFourthOrderTensorValue(data::AbstractArray{T}) where {T} = ((D,)=size(data); SymFourthOrderTensorValue{D}(data))
SymFourthOrderTensorValue{D}(data::AbstractArray{T}) where {D,T} = SymFourthOrderTensorValue{D,T}(_flatten_sym_fourth_order_tensor(data,Val(D)))
SymFourthOrderTensorValue{D,T1}(data::AbstractArray{T2}) where {D,T1,T2} = SymFourthOrderTensorValue{D,T1}(_flatten_sym_fourth_order_tensor(data,Val(D)))
SymFourthOrderTensorValue{D,T1,L}(data::AbstractArray{T2}) where {D,T1,T2,L} = SymFourthOrderTensorValue{D,T1,L}(_flatten_sym_fourth_order_tensor(data,Val(D)))

###############################################################
# Conversions (SymFourthOrderTensorValue)
###############################################################

@generated function _SymFourthOrder_to_array(arg::SymFourthOrderTensorValue{D,T,L}) where {D,T,L}
  str = ""
  for l in 1:D
    for k in 1:D
      for j in 1:D
        for i in 1:D
          str *= "arg[$i,$j,$k,$l], "
        end
      end
    end
  end
  Meta.parse("SArray{Tuple{D,D,D,D},T}(($str))")
end

# Inverse conversion
convert(::Type{<:MArray{Tuple{D,D,D,D},T}}, arg::SymFourthOrderTensorValue) where {D,T} = MArray{Tuple{D,D,D,D},T}(_SymFourthOrder_to_array(arg))
convert(::Type{<:SArray{Tuple{D,D,D,D},T}}, arg::SymFourthOrderTensorValue) where {D,T} = _SymFourthOrder_to_array(arg)

# Internal conversion
convert(::Type{<:SymFourthOrderTensorValue{D,T}}, arg::SymFourthOrderTensorValue{D}) where {D,T} = SymFourthOrderTensorValue{D,T}(Tuple(arg))
convert(::Type{<:SymFourthOrderTensorValue{D,T}}, arg::SymFourthOrderTensorValue{D,T}) where {D,T} = arg

###############################################################
# Other constructors and conversions (SymFourthOrderTensorValue)
###############################################################

# This is in fact the "symmetrized" 4th order identity
"""
    one(::SymFourthOrderTensorValue{D,T}})

Returns the tensor `resᵢⱼₖₗ = δᵢₖδⱼₗ(δᵢⱼ + (1-δᵢⱼ)/2)`.

The scalar type `T2` of the result is `typeof(one(T)/2)`.
"""
@generated function one(::Type{<:SymFourthOrderTensorValue{D,T}}) where {D,T}
  S = typeof(one(T)/2)
  str = join(["($i==$k && $j==$l) ?  ( $i==$j ? one($S) :  one($S)/2) : zero($S), " for i in 1:D for j in i:D for k in 1:D for l in k:D])
  Meta.parse("SymFourthOrderTensorValue{D,$S}(($str))")
end

change_eltype(::Type{<:SymFourthOrderTensorValue{D}},::Type{T2}) where {D,T2} = SymFourthOrderTensorValue{D,T2}
change_eltype(::Type{SymFourthOrderTensorValue{D,T1,L}},::Type{T2}) where {D,T1,T2,L} = SymFourthOrderTensorValue{D,T2,L}

num_indep_components(::Type{<:SymFourthOrderTensorValue{D}}) where {D} = (D*(D+1)÷2)^2

###############################################################
# VTK export (SymFourthOrderTensorValue)
###############################################################

function indep_components_names(::Type{<:SymFourthOrderTensorValue{A}}) where A
  if A>3
    return ["$i$j$k$l" for i in 1:A for j in i:A for k in 1:A for l in k:A ]
  else
    c_name = ["X", "Y", "Z"]
    return [c_name[i]*c_name[j]*c_name[k]*c_name[l] for i in 1:A for j in i:A for k in 1:A for l in k:A ]
  end
end

