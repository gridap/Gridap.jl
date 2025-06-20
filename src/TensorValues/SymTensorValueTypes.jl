###############################################################
# SymTensorValue Type
###############################################################
"""
    AbstractSymTensorValue{D,T,L} <: MultiValue{Tuple{D,D},T,2,L}

Abstract type representing any symmetric second-order `D`×`D` tensor, with symmetry ij↔ji.

See also [`SymTensorValue`](@ref), [`SymTracelessTensorValue`](@ref).
"""
abstract type AbstractSymTensorValue{D,T,L} <: MultiValue{Tuple{D,D},T,2,L} end

num_indep_components(::Type{<:AbstractSymTensorValue}) = @unreachable "Concrete
  type and dimension are needed to count components of symmetric tensors."

"""
    SymTensorValue{D,T,L} <: AbstractSymTensorValue{D,T,L}

Type representing a symmetric second-order `D`×`D` tensor. It must hold `L` = `D`(`D`+1)/2.

It is constructed by providing the components of index (i,j) for 1 ≤ i ≤ j ≤ `D`.
"""
struct SymTensorValue{D,T,L} <: AbstractSymTensorValue{D,T,L}
    data::NTuple{L,T}
    function SymTensorValue{D,T}(data::NTuple{L,T}) where {D,T,L}
        @check L == D*(D+1)÷2
        new{D,T,L}(data)
    end
end

function promote_rule(::Type{<:SymTensorValue{D,Ta}}, ::Type{<:SymTensorValue{D,Tb}}) where {D,Ta,Tb}
    T = promote_type(Ta,Tb)
    SymTensorValue{D,T}
end

###############################################################
# Constructors (SymTensorValue)
###############################################################

# Empty SymTensorValue constructor

SymTensorValue()                   = SymTensorValue{0,Int}(NTuple{0,Int}())
SymTensorValue{0}()                = SymTensorValue{0,Int}(NTuple{0,Int}())
SymTensorValue{0,T}() where {T}    = SymTensorValue{0,T}(NTuple{0,T}())
SymTensorValue(data::NTuple{0})    = SymTensorValue{0,Int}(data)
SymTensorValue{0}(data::NTuple{0}) = SymTensorValue{0,Int}(data)

# SymTensorValue single NTuple argument constructor

@generated function SymTensorValue(data::NTuple{L,T}) where {L,T}
  msg = "Invalid number of scalar arguments in SymTensorValue constructor"
  V = (sqrt(1+8*L)-1)/2
  @check floor(Int,V) == ceil(Int,V) msg
  D = Int(V)
  quote
    SymTensorValue{$D,T}(data)
  end
end
SymTensorValue{D}(data::NTuple{L,T}) where {D,L,T}           = SymTensorValue{D,T}(data)
SymTensorValue{D,T1}(data::NTuple{L,T2}) where {D,L,T1,T2}   = SymTensorValue{D,T1}(NTuple{L,T1}(data))
SymTensorValue{D,T1,L}(data::NTuple{L,T2}) where {D,L,T1,T2} = SymTensorValue{D,T1}(NTuple{L,T1}(data))

# SymTensorValue single Tuple argument constructor

SymTensorValue(data::Tuple)                      = SymTensorValue(promote(data...))
SymTensorValue{D}(data::Tuple) where {D}         = SymTensorValue{D}(promote(data...))
SymTensorValue{D,T1}(data::Tuple) where {D,T1}   = SymTensorValue{D,T1}(NTuple{length(data),T1}(data))
SymTensorValue{D,T1,L}(data::Tuple) where {D,T1,L}   = SymTensorValue{D,T1}(NTuple{L,T1}(data))

# SymTensorValue Vararg constructor

SymTensorValue(data::Number...)                    = SymTensorValue(data)
SymTensorValue{D}(data::Number...)    where {D}    = SymTensorValue{D}(data)
SymTensorValue{D,T1}(data::Number...) where {D,T1} = SymTensorValue{D,T1}(data)
SymTensorValue{D,T1,L}(data::Number...) where {D,T1,L} = SymTensorValue{D,T1}(data)

# SymTensorValue single AbstractMatrix argument constructor

#From Square Matrices
@generated function _flatten_upper_triangle(data::AbstractArray,::Val{D}) where D
  check_e = :( @check ($D,$D) == size(data) )
  str = ""
  for i in 1:D
    for j in i:D
      str *= "data[$i,$j], "
    end
  end
  ret_e = Meta.parse(" return ($str)")
  Expr(:block, check_e, ret_e)
end

SymTensorValue(data::AbstractMatrix{T}) where {T} = ((D1,)=size(data);  SymTensorValue{D1}(data))
SymTensorValue{D}(data::AbstractMatrix{T}) where {D,T} = SymTensorValue{D,T}(_flatten_upper_triangle(data,Val{D}()))
SymTensorValue{D,T1}(data::AbstractMatrix{T2}) where {D,T1,T2} = SymTensorValue{D,T1}(_flatten_upper_triangle(data,Val{D}()))
SymTensorValue{D,T1,L}(data::AbstractMatrix{T2}) where {D,T1,T2,L} = SymTensorValue{D,T1,L}(_flatten_upper_triangle(data,Val{D}()))

###############################################################
# Conversions (SymTensorValue)
###############################################################

@generated function _SymTensorValue_to_array(arg::SymTensorValue{D,T,L}) where {D,T,L}
  str = ""
  for j in 1:D
    for i in 1:D
      p = _2d_sym_tensor_linear_index(D,i,j)
      str *= "arg.data[$p], "
    end
  end
  Meta.parse("SMatrix{D,D,T}(($str))")
end

# Inverse conversion
convert(::Type{<:MArray{Tuple{D,D},T}}, arg::SymTensorValue) where {D,T} = MMatrix{D,D,T}(_SymTensorValue_to_array(arg))
convert(::Type{<:SArray{Tuple{D,D},T}}, arg::SymTensorValue) where {D,T} = _SymTensorValue_to_array(arg)

# Internal conversion
convert(::Type{<:SymTensorValue{D,T}}, arg::SymTensorValue{D}) where {D,T} = SymTensorValue{D,T}(Tuple(arg))
convert(::Type{<:SymTensorValue{D,T}}, arg::SymTensorValue{D,T}) where {D,T} = arg

###############################################################
# Other constructors and conversions (SymTensorValue)
###############################################################

@generated function one(::Type{<:SymTensorValue{D,T}}) where {D,T}
  str = join(["$i==$j ? one(T) : zero(T), " for i in 1:D for j in i:D])
  Meta.parse("SymTensorValue{D,T}(($str))")
end

change_eltype(::Type{<:SymTensorValue{D}},::Type{T2}) where {D,T2} = SymTensorValue{D,T2}
change_eltype(::Type{SymTensorValue{D,T1,L}},::Type{T2}) where {D,T1,T2,L} = SymTensorValue{D,T2,L}

num_indep_components(::Type{<:SymTensorValue{D}}) where {D} = D*(D+1)÷2

###############################################################
# VTK export (SymTensorValue)
###############################################################

function indep_components_names(::Type{<:AbstractSymTensorValue{D}}) where D
  if D>3
    return ["$i$j" for i in 1:D for j in i:D ]
  else
    c_name = ["X", "Y", "Z"]
    return [c_name[i]*c_name[j] for i in 1:D for j in i:D ]
  end
end

