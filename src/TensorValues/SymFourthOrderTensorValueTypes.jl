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

###############################################################
# Conversions (SymFourthOrderTensorValue)
###############################################################

# Direct conversion
convert(::Type{<:SymFourthOrderTensorValue{D,T}}, arg::Tuple) where {D,T} = SymFourthOrderTensorValue{D,T}(arg)

# Inverse conversion
convert(::Type{<:NTuple{L,T}}, arg::SymFourthOrderTensorValue) where {L,T} = NTuple{L,T}(Tuple(arg))

# Internal conversion
convert(::Type{<:SymFourthOrderTensorValue{D,T}}, arg::SymFourthOrderTensorValue{D}) where {D,T} = SymFourthOrderTensorValue{D,T}(Tuple(arg))
convert(::Type{<:SymFourthOrderTensorValue{D,T}}, arg::SymFourthOrderTensorValue{D,T}) where {D,T} = arg

###############################################################
# Other constructors and conversions (SymFourthOrderTensorValue)
###############################################################

@generated function zero(::Type{<:SymFourthOrderTensorValue{D,T}}) where {D,T}
  L=(D*(D+1)÷2)^2
  quote
    SymFourthOrderTensorValue{D,T}(tfill(zero(T),Val{$L}()))
  end
end
zero(::Type{<:SymFourthOrderTensorValue{D,T,L}}) where {D,T,L} = SymFourthOrderTensorValue{D,T}(tfill(zero(T),Val{L}()))
zero(::SymFourthOrderTensorValue{D,T,L}) where {D,T,L} = zero(SymFourthOrderTensorValue{D,T,L})

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
one(::SymFourthOrderTensorValue{D,T}) where {D,T} = one(SymFourthOrderTensorValue{D,T})

@generated function rand(rng::AbstractRNG,
                         ::Random.SamplerType{<:SymFourthOrderTensorValue{D,T}}) where {D,T}
  L=(D*(D+1)÷2)^2
  quote
    rand(rng, SymFourthOrderTensorValue{D,T,$L})
  end
end
rand(rng::AbstractRNG,::Random.SamplerType{<:SymFourthOrderTensorValue{D,T,L}}) where {D,T,L} =
  SymFourthOrderTensorValue{D,T}(Tuple(rand(rng, SVector{L,T})))

Mutable(::Type{<:SymFourthOrderTensorValue{D,T}}) where {D,T} = @notimplemented
Mutable(::SymFourthOrderTensorValue{D,T}) where {D,T} = Mutable(SymFourthOrderTensorValue{D,T})
mutable(a::SymFourthOrderTensorValue{D}) where D = @notimplemented

change_eltype(::Type{SymFourthOrderTensorValue{D,T1}},::Type{T2}) where {D,T1,T2} = SymFourthOrderTensorValue{D,T2}
change_eltype(::Type{SymFourthOrderTensorValue{D,T1,L}},::Type{T2}) where {D,T1,T2,L} = SymFourthOrderTensorValue{D,T2,L}
change_eltype(::SymFourthOrderTensorValue{D,T1,L},::Type{T2}) where {D,T1,T2,L} = change_eltype(SymFourthOrderTensorValue{D,T1,L},T2)

###############################################################
# Introspection (SymFourthOrderTensorValue)
###############################################################

eltype(::Type{<:SymFourthOrderTensorValue{D,T}}) where {D,T} = T
eltype(::SymFourthOrderTensorValue{D,T}) where {D,T} = eltype(SymFourthOrderTensorValue{D,T})

size(::Type{<:SymFourthOrderTensorValue{D}}) where {D} = (D,D,D,D)
size(::SymFourthOrderTensorValue{D}) where {D} = size(SymFourthOrderTensorValue{D})

length(::Type{<:SymFourthOrderTensorValue{D}}) where {D} = D*D*D*D
length(::SymFourthOrderTensorValue{D}) where {D} = length(SymFourthOrderTensorValue{D})

num_components(::Type{<:SymFourthOrderTensorValue}) = @unreachable "The dimension is needed to count components"
num_components(::Type{<:SymFourthOrderTensorValue{D}}) where {D} = length(SymFourthOrderTensorValue{D})
num_components(::SymFourthOrderTensorValue{D}) where {D} = num_components(SymFourthOrderTensorValue{D})

num_indep_components(::Type{<:SymFourthOrderTensorValue})  = num_components(SymFourthOrderTensorValue)
num_indep_components(::Type{<:SymFourthOrderTensorValue{D}}) where {D} = (D*(D+1)÷2)^2
num_indep_components(::SymFourthOrderTensorValue{D}) where {D} = num_indep_components(SymFourthOrderTensorValue{D})

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

