###############################################################
# SymTracelessTensorValue Type
###############################################################

"""
    SymTracelessTensorValue{D,T,L} <: AbstractSymTensorValue{D,T,L}
    QTensorValue{D,T,L}

Type representing a symmetric second-order `D`×`D` tensor with zero trace. It must hold `L` = `D`(`D`+1)/2.
This type is used to model the Q-tensor order parameter in nematic liquid cristals.

The constructor determines the value of index (`D`,`D`) as minus the sum of the other diagonal values, so it value musn't be provided. The constructor thus expects the `L`-1 components of indices (i,j) for 1 ≤ i ≤ `D`-1 and  i ≤ j ≤ `D`.
"""
struct SymTracelessTensorValue{D,T,L} <: AbstractSymTensorValue{D,T,L}
  data::NTuple{L,T}
  function SymTracelessTensorValue{D,T}(data::NTuple{L,T}) where {D,T,L}
    @check L == num_indep_components(SymTracelessTensorValue{D})
    Int(D) == 0 && return new{D,T,0}(data)
    Int(D) == 1 && return new{D,T,1}( (zero(T),) )
    new{D,T,L+1}( (data..., _minus_trace(data,Val(D))) )
  end
end

@generated function _minus_trace(data::NTuple{L,T},::Val{D}) where {D,T,L}
  trace_expr = Expr(:call, :+)
  for i in 1:D-1
    k = _2d_sym_tensor_linear_index(D,i,i)
    push!(trace_expr.args, :(- data[$k]) )
  end
  return trace_expr
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

SymTracelessTensorValue()                = SymTracelessTensorValue{0,Int}(NTuple{0,Int}())
SymTracelessTensorValue{0,T}() where {T} = SymTracelessTensorValue{0,T}(NTuple{0,T}())
SymTracelessTensorValue(data::NTuple{0}) = SymTracelessTensorValue{0,Int}(data)

# A traceless tensor of dimension 0 or 1 is zero, so no component is given, for
# any variance of its dimension

SymTracelessTensorValue{D}()                where D = SymTracelessTensorValue{D,Int}(NTuple{0,Int}())
SymTracelessTensorValue{D}(data::NTuple{0}) where D = SymTracelessTensorValue{D,Int}(data)

# SymTracelessTensorValue single NTuple argument constructor

@generated function SymTracelessTensorValue(data::NTuple{L,T}) where {L,T}
  msg = "Invalid number of scalar arguments in SymTracelessTensorValue constructor"
  D = (isqrt(9+8*L)-1)÷2
  @check L == D*(D+1)÷2-1 msg
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
  comps = Expr[]
  for i in 1:D-1
    for j in i:D
      push!(comps, :( data[$i,$j]))
    end
  end

  quote
    @check ($D,$D) == size(data)
    return tuple($(comps...))
  end
end

SymTracelessTensorValue(data::AbstractMatrix{T}) where {T} = ((D1,D2)=size(data); SymTracelessTensorValue{D1}(data))
SymTracelessTensorValue{D}(data::AbstractMatrix{T}) where {D,T} = SymTracelessTensorValue{D,T}(_flatten_upper_triangle_traceless(data,Val{D}()))
SymTracelessTensorValue{D,T1}(data::AbstractMatrix{T2}) where {D,T1,T2} = SymTracelessTensorValue{D,T1}(_flatten_upper_triangle_traceless(data,Val{D}()))
SymTracelessTensorValue{D,T1,L}(data::AbstractMatrix{T2}) where {D,T1,T2,L} = SymTracelessTensorValue{D,T1,L}(_flatten_upper_triangle_traceless(data,Val{D}()))

###############################################################
# Conversions (SymTracelessTensorValue)
###############################################################

@generated function _SymTracelessTensorValue_to_array(arg::SymTracelessTensorValue{Da,T,L}) where {Da,T,L}
  D = Int(Da)
  comps = Expr[]
  for j in 1:D
    for i in 1:D
      push!(comps, :( arg[$i,$j]))
    end
  end
  :( SMatrix{$D,$D,T}($(comps...)) )
end

# Inverse conversion
convert(::Type{<:MArray{Tuple{D,D},T}}, arg::SymTracelessTensorValue) where {D,T} = MMatrix{D,D,T}(_SymTracelessTensorValue_to_array(arg))
convert(::Type{<:SArray{Tuple{D,D},T}}, arg::SymTracelessTensorValue) where {D,T} = SMatrix{D,D,T}(_SymTracelessTensorValue_to_array(arg))

# Internal conversion
convert(::Type{<:SymTracelessTensorValue{D}}, arg::SymTracelessTensorValue{D}) where {D} = SymTracelessTensorValue{D}(Tuple(arg)[1:end-1])
convert(::Type{<:SymTracelessTensorValue{D,T}}, arg::SymTracelessTensorValue{D}) where {D,T} = SymTracelessTensorValue{D,T}(Tuple(arg)[1:end-1])

###############################################################
# Other constructors and conversions (SymTracelessTensorValue)
###############################################################

change_eltype(::Type{<:SymTracelessTensorValue{D}},::Type{T2}) where {D,T2} = SymTracelessTensorValue{D,T2}
change_eltype(::Type{SymTracelessTensorValue{D,T1,L}},::Type{T2}) where {D,T1,T2,L} = SymTracelessTensorValue{D,T2,L}

function change_variance(::Type{<:SymTracelessTensorValue{D}}, ::Val{Sv}) where {D,Sv}
  _check_change_variance(SymTracelessTensorValue{D}, Sv, true)
  SymTracelessTensorValue{first(Sv)}
end
function change_variance(::Type{<:SymTracelessTensorValue{D,T}}, ::Val{Sv}) where {D,T,Sv}
  _check_change_variance(SymTracelessTensorValue{D}, Sv, true)
  SymTracelessTensorValue{first(Sv),T}
end

num_indep_components(::Type{<:SymTracelessTensorValue{D}}) where {D} =
  iszero(Int(D)) ? 0 : D*(D+1)÷2-1

num_stored_components(::Type{<:SymTracelessTensorValue{D}}) where D = num_indep_components(SymTensorValue{D})

get_indep_components(a::SymTracelessTensorValue) = Base.front(Tuple(a))

stored_comp_indices(::Type{<:SymTracelessTensorValue{D}}) where D = stored_comp_indices(SymTensorValue{D})

