# ExteriorFormValues.jl

"""
    SkewSymKTensorValue{V,K,D,T,L} <: MultiValue{NTuple{K,D},T,K,L}

Value of a completely skew-symmetric (alternating) `K`-th order tensor in `D`
dimensions. It stores `L = binomial(D,K)` independent components, indexed by
the K-combinations `I` of `1:D`, in lexicographic order.

All `K` indices share the variance `V<:Dim`, which the constructor
automatically gives to `D`: `SkewSymKTensorValue{1,2,Int,2,Cov}(1,2)` works but
returns a `SkewSymKTensorValue{1,cov2,Int,2,Cov}`.

Use the aliases [`ExteriorFormValue`](@ref) and [`ExteriorVectorValue`](@ref)
rather than this type; they fix `V` and are the types the constructors,
`convert` and `show` are defined on.
"""
struct SkewSymKTensorValue{V,K,D,T,L} <: MultiValue{NTuple{K,D},T,K,L}
  data::NTuple{L,T}

  function SkewSymKTensorValue{V,K,D,T,L}(data::Tuple) where {V<:Dim,K,D,T,L}
    @assert D isa Int || D isa V  "implicit conversion of $D::$(typeof(D)) to $V is not allowed"
    @assert L == binomial(Int(D),K) "wrong number of independent components: got $L, expected $(binomial(Int(D),K))"
    @assert length(data) == L "wrong number of values: got $(length(data)), expected $L"
    new{V,K,V(Int(D)),T,L}(NTuple{L,T}(data))
  end
end

"""
    ExteriorFormValue{K,D,T,L} = SkewSymKTensorValue{Cov,K,D,T,L}

Alias for an alternating `K`-form in `D` dimensions with components in the
basis `{dx^I}`.

`show` labels the components `dxⁱ`, or `dλⁱ` when the `IOContext` property
`:coordinates` is `:barycentric`.
"""
const ExteriorFormValue{K,D,T,L} = SkewSymKTensorValue{Cov,K,D,T,L}

"""
    ExteriorVectorValue{K,D,T,L} = SkewSymKTensorValue{Con,K,D,T,L}

Alias for an alternating `K`-vector in `D` dimensions with components in the
basis `{∂_I}`.

`show` labels the components `∂ᵢ`, or `∂λᵢ` when the `IOContext` property
`:coordinates` is `:barycentric`.
"""
const ExteriorVectorValue{K,D,T,L} = SkewSymKTensorValue{Con,K,D,T,L}
#const _UnvariancedSkewSymKTensorValue{K,D,T,L} = SkewSymKTensorValue{Int,K,D,T,L}

###############################################################
# Constructors (SkewSymKTensorValue)
###############################################################

SkewSymKTensorValue{V,K,D,T,L}(data::Number...) where {V,K,D,T,L} = SkewSymKTensorValue{V,K,D,T,L}(data)
SkewSymKTensorValue{V,K,D,T}(data::NTuple{L,T}) where {V,K,D,T,L} = SkewSymKTensorValue{V,K,D,T,L}(data)

# Empty constructors
SkewSymKTensorValue{V,K,D}()                where {V,K,D}   = SkewSymKTensorValue{V,K,D,Int}(NTuple{0,Int}())
SkewSymKTensorValue{V,K,D,T}()              where {V,K,D,T} = SkewSymKTensorValue{V,K,D,T}(NTuple{0,T}())
SkewSymKTensorValue{V,K,D}(data::NTuple{0}) where {V,K,D}   = SkewSymKTensorValue{V,K,D,Int}(data)

# NTuple argument constructors
SkewSymKTensorValue{V,K,D}(data::NTuple{L,T}) where {V,K,D,L,T}         = SkewSymKTensorValue{V,K,D,T}(data)
SkewSymKTensorValue{V,K,D,T1}(data::NTuple{L,T2}) where {V,K,D,L,T1,T2} = SkewSymKTensorValue{V,K,D,T1}(NTuple{L,T1}(data))

# single Tuple argument constructors
SkewSymKTensorValue{V,K,D}(data::Tuple) where {V,K,D}       = SkewSymKTensorValue{V,K,D}(promote(data...))
SkewSymKTensorValue{V,K,D,T1}(data::Tuple) where {V,K,D,T1} = SkewSymKTensorValue{V,K,D,T1}(NTuple{length(data),T1}(data))

# Vararg constructors
SkewSymKTensorValue{V,K,D}(data::Number...) where {V,K,D}       = SkewSymKTensorValue{V,K,D}(data)
SkewSymKTensorValue{V,K,D,T1}(data::Number...) where {V,K,D,T1} = SkewSymKTensorValue{V,K,D,T1}(data)

######################################
# Conversions and other constructors #
######################################

function convert(::Type{<:SkewSymKTensorValue{V,K,D,T}}, arg::SkewSymKTensorValue{V,K,Da}) where {V,K,D,Da,T}
  D == Da || error("Cannot convert value of type $(typeof(arg)) to type SkewSymKTensorValue{$V,$K,$D,$T}")
  SkewSymKTensorValue{V,K,D,T}(arg.data)
end
function convert(::Type{<:SkewSymKTensorValue{V,K,D,T}}, arg::SkewSymKTensorValue{V,K,Da,T}) where {V,K,D,Da,T}
  D == Da || error("Cannot convert value of type $(typeof(arg)) to type SkewSymKTensorValue{$V,$K,$D,$T}")
  SkewSymKTensorValue{V,K,D,T}(arg.data)
end

for (Alias,V) in ((:ExteriorFormValue,:Cov), (:ExteriorVectorValue,:Con))
  @eval begin
    # Internal conversion, which never changes the variance
    function convert(::Type{<:$Alias{K,D,T}}, arg::$Alias{K,Da}) where {K,D,Da,T}
      D == Da || error("Cannot convert value of type $(typeof(arg)) to type $($Alias){$K,$D,$T}")
      SkewSymKTensorValue{V,K,D,T}(arg.data)
    end
    function convert(::Type{<:$Alias{K,D,T}}, arg::$Alias{K,Da,T}) where {K,D,Da,T}
      D == Da || error("Cannot convert value of type $(typeof(arg)) to type $($Alias){$K,$D,$T}")
      SkewSymKTensorValue{V,K,D,T}(arg.data)
    end
  end
end

# Full tensor components:
# - `A[I] = 0` when `I` repeats an index, else
# - `A[I] = levicivita(σ) * ω_{sort(I)}`, where `σ` sorts `I`
@generated function _SkewSymKTensorValue_to_array(arg::SkewSymKTensorValue{V,K,Da,T,L}) where {V,K,Da,T,L}
  D = Int(Da)
  comps = Expr[]
  for idx in CartesianIndices(ntuple(_ -> D, Val(K)))
    I = Tuple(idx)
    s = sorting_sign(I...)
    if s == 0
      push!(comps, :(zero(T)))
    else
      l = combination_index(sort(SVector{K,Int}(I)), D)
      push!(comps, s > 0 ? :(arg.data[$l]) : :(-arg.data[$l]))
    end
  end

  S = Tuple{ntuple(_ -> D, K)...}
  quote
    return SArray{$S,T,K,$(D^K)}( tuple($(comps...)) )
  end
end

# Inverse conversion
convert(::Type{<:MArray{S,T}}, arg::SkewSymKTensorValue) where {S,T} = MArray{S,T}(_SkewSymKTensorValue_to_array(arg))
convert(::Type{<:SArray{S,T}}, arg::SkewSymKTensorValue) where {S,T} = SArray{S,T}(_SkewSymKTensorValue_to_array(arg))

one(::Type{<:SkewSymKTensorValue{V,0,D,T}}) where {D,T} = one(T)
one(::Type{<:SkewSymKTensorValue{V,K}}) where K = @unreachable "Skew symmetric k-tensors do not have multiplicative neutral for k > 0."

change_eltype(::Type{<:SkewSymKTensorValue{V,K,D,T1,L}}, ::Type{T2}) where {V,K,D,T1,T2,L} =
  SkewSymKTensorValue{V,K,D,T2,L}

num_indep_components(::Type{<:SkewSymKTensorValue{V,K,D}}) where {V,K,D} = binomial(Int(D),K)

# ============================================================
# Display
# ============================================================

const _sbs_cart_cov = ["dx¹","dx²","dx³","dx⁴","dx⁵","dx⁶","dx⁷","dx⁸","dx⁹"]
const _sbs_bary_cov = ["dλ¹","dλ²","dλ³","dλ⁴","dλ⁵","dλ⁶","dλ⁷","dλ⁸","dλ⁹"]
const _sbs_cart_con = ["∂₁","∂₂","∂₃","∂₄","∂₅","∂₆","∂₇","∂₈","∂₉"]
const _sbs_bary_con = ["∂λ₁","∂λ₂","∂λ₃","∂λ₄","∂λ₅","∂λ₆","∂λ₇","∂λ₈","∂λ₉"]

# Labels of the frame the components of a `V`-variant tensor sit on: the
# coframe `{dxⁱ}` for forms, the frame `{∂ᵢ}` dual to it for vectors.
_frame_labels(::Type{Cov}, ::Val{:cartesian})   = _sbs_cart_cov
_frame_labels(::Type{Cov}, ::Val{:barycentric}) = _sbs_bary_cov
_frame_labels(::Type{Con}, ::Val{:cartesian})   = _sbs_cart_con
_frame_labels(::Type{Con}, ::Val{:barycentric}) = _sbs_bary_con

function _frame_labels(io::IO, ::Type{V}) where V
  c = get(io, :coordinates, :cartesian)
  (c === :cartesian || c === :barycentric) || throw(ArgumentError(
    "unknown :coordinates value $(repr(c)); expected :cartesian or :barycentric"))
  _frame_labels(V, Val(c))
end

# `V` is not deducible from a partially parametrised type, so the frame is
# picked per alias, where it is fixed.
indep_components_names(::Type{<:ExteriorFormValue{K,D}}) where {K,D} =
  [join(_frame_labels(Cov, Val(:cartesian))[I], " ∧ ") for I in sorted_combinations(Int(D), K)]
indep_components_names(::Type{<:ExteriorVectorValue{K,D}}) where {K,D} =
  [join(_frame_labels(Con, Val(:cartesian))[I], " ∧ ") for I in sorted_combinations(Int(D), K)]

function _show_skewsymk(io::IO, a::SkewSymKTensorValue{V,K,D}, basis_strs) where {V,K,D}
  L = length(a.data)
  L == 0 && return print(io, "0")
  cs = sorted_combinations(Int(D), K)
  parts = [join(basis_strs[I], " ∧ ") for I in cs]
  result = [string(ai, " ", bi) for (ai, bi) in zip(a.data, parts)]
  print(io, join(result, " + "))
end

"""
    show(io::IO, ::MIME"text/plain", ω::SkewSymKTensorValue)

Print `ω` as a linear combination of frame elements, labelled `dx¹,…,dx^D` for
an [`ExteriorFormValue`](@ref) and `∂₁,…,∂_D` for an
[`ExteriorVectorValue`](@ref).

Set the `:coordinates` `IOContext` property to `:barycentric` to label them
`dλ¹,…,dλ^D` (resp. `∂λ₁,…,∂λ_D`) instead, for a tensor whose components are
coefficients on the ambient barycentric (co)frame of a (D−1)-simplex:

    show(IOContext(stdout, :coordinates => :barycentric), MIME("text/plain"), ω)
"""
function Base.show(io::IO, ::MIME"text/plain", a::SkewSymKTensorValue{V,K,D,T,L}) where {V,K,D,T,L}
  @assert D <= 9 "show not implemented for D > 9"
  _show_skewsymk(io, a, _frame_labels(io, V))
end

Base.show(io::IO, mime::MIME"text/plain", a::SkewSymKTensorValue{K,D,T,L,Int}) where {K,D,T,L} =
  @invoke show(io::IO, mime::MIME"text/plain", a::MultiValue)

