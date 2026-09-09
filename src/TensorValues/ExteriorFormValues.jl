# ExteriorFormValues.jl
#

"""
    ExteriorFormValue{K,D,T,L} <: MultiValue{NTuple{K,D},T,K,L}

Value of a differential K-form in D dimensions storing the `L = binomial(D,K)`
independent components on the orientation-ordered basis `{dx^I}` (lexicographic
K-combinations `I` of `1:D`), with scalar type `T`.

A form takes covariant indices, so `D` is stored as a [`Cov`](@ref) dimension:
`ExteriorFormValue{1,3}` and `ExteriorFormValue{1,cov3}` name the same type, and
`size` returns `(cov3,)`. Use `Int(D)` wherever the dimension is needed as a
number.

`show` labels the components `dxⁱ`, or `dλⁱ` when the `IOContext` property
`:coordinates` is `:barycentric`.
"""
struct ExteriorFormValue{K,D,T,L} <: MultiValue{NTuple{K,D},T,K,L}
  data::NTuple{L,T}

  function ExteriorFormValue{K,D,T}(data::NTuple{L,T}) where {K,D,T,L}
    @assert L == binomial(Int(D),K) "wrong number of values: got $L, expected $(binomial(Int(D),K))"
    new{K,Cov(D),T,L}(data)
  end
end

###############################################################
# Constructors (ExteriorFormValue)
###############################################################

# Empty constructors
ExteriorFormValue{K,D}()                where {K,D}   = ExteriorFormValue{K,D,Int}(NTuple{0,Int}())
ExteriorFormValue{K,D,T}()              where {K,D,T} = ExteriorFormValue{K,D,T}(NTuple{0,T}())
ExteriorFormValue{K,D}(data::NTuple{0}) where {K,D}   = ExteriorFormValue{K,D,Int}(data)

# NTuple argument constructors
ExteriorFormValue{K,D}(data::NTuple{L,T}) where {K,D,L,T}           = ExteriorFormValue{K,D,T}(data)
ExteriorFormValue{K,D,T1}(data::NTuple{L,T2}) where {K,D,L,T1,T2}   = ExteriorFormValue{K,D,T1}(NTuple{L,T1}(data))
ExteriorFormValue{K,D,T1,L}(data::NTuple{L,T2}) where {K,D,L,T1,T2} = ExteriorFormValue{K,D,T1}(NTuple{L,T1}(data))

# single Tuple argument constructors
ExteriorFormValue{K,D}(data::Tuple) where {K,D}            = ExteriorFormValue{K,D}(promote(data...))
ExteriorFormValue{K,D,T1}(data::Tuple) where {K,D,T1}      = ExteriorFormValue{K,D,T1}(NTuple{length(data),T1}(data))
ExteriorFormValue{K,D,T1,L}(data::Tuple) where {K,D,T1,L}  = ExteriorFormValue{K,D,T1}(NTuple{L,T1}(data))

# Vararg constructors
ExteriorFormValue{K,D}(data::Number...) where {K,D}           = ExteriorFormValue{K,D}(data)
ExteriorFormValue{K,D,T1}(data::Number...) where {K,D,T1}     = ExteriorFormValue{K,D,T1}(data)
ExteriorFormValue{K,D,T1,L}(data::Number...) where {K,D,T1,L} = ExteriorFormValue{K,D,T1}(data)


######################################
# Conversions and other constructors #
######################################

# Full tensor components:
# - `A[I] = 0` when `I` repeats an index, else
# - `A[I] = levicivita(σ) * ω_{sort(I)}`, where `σ` sorts `I`
@generated function _FormValue_to_array(arg::ExteriorFormValue{K,Da,T,L}) where {K,Da,T,L}
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
convert(::Type{<:MArray{S,T}}, arg::ExteriorFormValue) where {S,T} = MArray{S,T}(_FormValue_to_array(arg))
convert(::Type{<:SArray{S,T}}, arg::ExteriorFormValue) where {S,T} = SArray{S,T}(_FormValue_to_array(arg))

# Internal conversion
function convert(::Type{<:ExteriorFormValue{K,D,T}}, arg::ExteriorFormValue{K,Da}) where {K,D,Da,T}
  D == Da || error("Cannot convert value of type $(typeof(arg)) to type ExteriorFormValue{$K,$D,$T}")
  ExteriorFormValue{K,D,T}(arg.data)
end
function convert(::Type{<:ExteriorFormValue{K,D,T}}, arg::ExteriorFormValue{K,Da,T}) where {K,D,Da,T}
  D == Da || error("Cannot convert value of type $(typeof(arg)) to type ExteriorFormValue{$K,$D,$T}")
  ExteriorFormValue{K,D,T}(arg.data)
end

one(::Type{<:ExteriorFormValue{0,D,T}}) where {D,T} = one(T)
one(::Type{<:ExteriorFormValue{K}}) where K = @unreachable "Exterior k-form do not have multiplicative neutral for k > 0."

change_eltype(::Type{<:ExteriorFormValue{K,D,T1,L}}, ::Type{T2}) where {K,D,T1,T2,L} =
  ExteriorFormValue{K,D,T2,L}

num_indep_components(::Type{<:ExteriorFormValue{K,D}}) where {K,D} = binomial(Int(D),K)

# ============================================================
# Display
# ============================================================

const _sbs_cart = ["dx¹","dx²","dx³","dx⁴","dx⁵","dx⁶","dx⁷","dx⁸","dx⁹"]
const _sbs_bary = ["dλ¹","dλ²","dλ³","dλ⁴","dλ⁵","dλ⁶","dλ⁷","dλ⁸","dλ⁹"]

function indep_components_names(::Type{<:ExteriorFormValue{K,D}}) where {K,D}
  [join(_sbs_cart[I], " ∧ ") for I in sorted_combinations(Int(D), K)]
end

function _show_dfv(io::IO, a::ExteriorFormValue{K,D}, basis_strs) where {K,D}
  L = length(a.data)
  L == 0 && return print(io, "0")
  cs = sorted_combinations(Int(D), K)
  parts = [join(basis_strs[I], " ∧ ") for I in cs]
  result = [string(ai, " ", bi) for (ai, bi) in zip(a.data, parts)]
  print(io, join(result, " + "))
end

function _coframe_labels(io::IO)
  c = get(io, :coordinates, :cartesian)
  c === :cartesian   && return _sbs_cart
  c === :barycentric && return _sbs_bary
  throw(ArgumentError(
    "unknown :coordinates value $(repr(c)); expected :cartesian or :barycentric"))
end

"""
    show(io::IO, ::MIME"text/plain", ω::ExteriorFormValue)

Print `ω` as a linear combination of coframe elements, labelled `dx¹,…,dx^D`.

Set the `:coordinates` `IOContext` property to `:barycentric` to label them
`dλ¹,…,dλ^D` instead, for a form whose components are coefficients on the
ambient barycentric coframe of a (D−1)-simplex:

    show(IOContext(stdout, :coordinates => :barycentric), MIME("text/plain"), ω)
"""
function Base.show(io::IO, ::MIME"text/plain", a::ExteriorFormValue{K,D}) where {K,D}
  @assert D <= 9 "show not implemented for D > 9"
  _show_dfv(io, a, _coframe_labels(io))
end

