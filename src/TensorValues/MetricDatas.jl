# Metric data handler

"""
    MetricData(::UniformScaling)
    MetricData(g, g_inv, sqrt_abs_g)

Helper struct to handle the metric information consisting of the metric tensor
`g`, its inverse `g_inv`, and `sqrt_abs_g`=√|det(g)|.

See [`parse_metric_data`](@ref) for lazy computation, and
[`compute_metric_data`](@ref) for greedy computation.

```
_g = SymTensorValue{cov2}(1.,2,3)
_g_inv = inv(_g)
_sqrt_abs_g = sqrt(abs(det(_g)))

g_or_ginv      = rand(Bool) ? _g : _g_inv
det_or_nothing = rand(Bool) ? _sqrt_abs_g : nothing

md = compute_metric_data(g_or_ginv, det_or_nothing)
md.g ≈ _g && md.g_inv ≈ _g_inv && md.sqrt_abs_g ≈ _sqrt_abs_g # true

md = parse_metric_data(g_or_ginv, det_or_nothing)
(; g, g_inv, sqrt_abs_g) = md # struct destructuring, lazily computing destructured fields only
   g ≈ _g &&    g_inv ≈ _g_inv &&    sqrt_abs_g ≈ _sqrt_abs_g  # true
```
"""
struct MetricData{
  # G and Ginv are both Nothing if and only if the metric is UniformScaling{T}(one(T))
  # non-nothing G and Ginv are covariant and contravariant, respectively
  G    <: Union{SymTensorValue,Nothing},
  Ginv <: Union{SymTensorValue,Nothing},
  T    <: Union{Number,Nothing}
}
  g::G
  g_inv::Ginv
  sqrt_abs_g::T

  function MetricData(I::UniformScaling{T}) where T
    @assert I.λ == one(T) "UniformScaling can only be used as a metric if it is the scaling by 1"
    new{Nothing,Nothing,T}(nothing,nothing,one(T))
  end

  function MetricData(g::Union{SymTensorValue{Dcov,Tg},Nothing},
                      g_inv::Union{SymTensorValue{Dcon,Ti},Nothing},
                      sq_g::Ts) where {Dcov,Tg,Dcon,Ti,Ts}

    isnothing(g) && isnothing(g_inv) && throw(ArgumentError("either the metric tensor g or its inverse g_inv must be provided, got neither"))

    # validate dimension and variance
    _Dcov = isnothing(g)     ? inverse_variance(Dcon) : Dcov
    _Dcon = isnothing(g_inv) ? inverse_variance(Dcov) : Dcon
    @assert _Dcov == _Dcon "metric tensor g and inverse g_inv do not have the same dimension"
    @assert _Dcov isa Cov  "metric tensor g must be covariant"
    @assert _Dcon isa Con  "inverse metric tensor g_inv must be contravariant"

    # promote to same scalar type
    _Tg = isnothing(g)     ? Bool : Tg
    _Ti = isnothing(g_inv) ? Bool : Ti
    _Ts = isnothing(sq_g)  ? Bool : Ts
    _T = Base.promote_op(inv, promote_type(_Tg, _Ti, _Ts))
    _g     = isnothing(g)     ? g     : convert(SymTensorValue{_Dcov,_T}, g)
    _g_inv = isnothing(g_inv) ? g_inv : convert(SymTensorValue{_Dcon,_T}, g_inv)
    _sq_g  = isnothing(sq_g)  ? sq_g  : convert(_T, sq_g)

    @check(isnothing(_g) || isnothing(_g_inv) || _g⋅_g_inv ≈ one(TensorValue{_Dcov,_Dcon,_T}),
           "the metric tensor g and its inverse g_inv are not inverse of each other")

    G, Ginv, T = typeof(_g), typeof(_g_inv), typeof(_sq_g)
    new{G,Ginv,T}(_g, _g_inv, _sq_g)
  end
end

# enable splating metric data in keyword arguments
Base.NamedTuple(md::MetricData) = (g=md.g, g_inv=md.g_inv, sqrt_abs_g=md.sqrt_abs_g)
Base.merge(nt::NamedTuple, md::MetricData) = merge(nt, NamedTuple(md))
MetricData(gdata::MetricData) = gdata

# Computes the field (only) if it is not stored
function Base.getproperty(md::MetricData, name::Symbol)
  if     name==:g
    g = getfield(md, :g)
    !isnothing(g) && return g
    return inv(getfield(md, :g_inv))

  elseif name==:g_inv
    g_inv = getfield(md, :g_inv)
    !isnothing(g_inv) && return g_inv
    return inv(getfield(md, :g))

  elseif name==:sqrt_abs_g
    sqrt_abs_g = getfield(md, :sqrt_abs_g)
    !isnothing(sqrt_abs_g) && return sqrt_abs_g
    g = getfield(md, :g)
    !isnothing(g) && return _sqrt_abs_det(g)
    g_inv = getfield(md, :g_inv)
    return inv(_sqrt_abs_det(g_inv))

  else
    throw(Base.FieldError(MetricData, name))
  end
end

function Base.getproperty(::MetricData{Nothing,Nothing,T}, name::Symbol) where T
  if     name==:g
    return UniformScaling{T}(one(T))
  elseif name==:g_inv
    return UniformScaling{T}(one(T))
  elseif name==:sqrt_abs_g
    return one(T)
  else
    throw(Base.FieldError(MetricData, name))
  end
end

# The volume density √|det g| of a metric tensor. A function of its own so
# that `partial_derivatives` of an expression holding it has a rule to apply,
# see [`partials_op_rule`](@ref).
_sqrt_abs_det(g::SymTensorValue) = sqrt(abs(det(g)))

function _check_metric_dim(g_or_ginv::SymTensorValue{Dg}, D::Int) where Dg
  Int(Dg) == D || throw(DimensionMismatch(
    "the metric tensor has dimension $(Int(Dg)), but the tensor it acts on has dimension $D"))
  nothing
end

"""
    isidentity(md::MetricData) -> Bool

Whether `md` represents the iddentity metric `LinearAlgebra.I`.
"""
isidentity(::MetricData{Nothing,Nothing}) = true
isidentity(::MetricData) = false

# Metric arguments: Type of valid slurped arguments to forward to
# parse_metric_data or compute_metric_data
const _MetricArgNoDetG = Union{MetricData,SymTensorValue,UniformScaling} # for norm(a, p, metric_args_nodet...)
const _MetricArg = Union{_MetricArgNoDetG,Number}


"""
    compute_metric_data(metric_args...)
    compute_metric_data(md::MetricData)
    compute_metric_data(LinearAlgebra.I)

Return a [`MetricData`](@ref) behaving like `(;g, g_inv, sqrt_abs_g)` with all
three fields eagerly computed. `metric_args` are the positional arguments
accepted by [`parse_metric_data`](@ref), in any order.
"""
compute_metric_data(metric_args...) = compute_metric_data(parse_metric_data(metric_args...))
compute_metric_data(I::UniformScaling) = MetricData(I)
compute_metric_data(md::MetricData{Nothing,Nothing}) = md # nothing to compute for the identity
compute_metric_data(md::MetricData) = MetricData(md.g, md.g_inv, md.sqrt_abs_g) # getproperty computes missing ones

# When a single tensor is stored, its inverse and √|det| share one determinant.
function compute_metric_data(md::MetricData{<:SymTensorValue,Nothing,Nothing})
  g = getfield(md, :g)
  det_g = det(g)
  MetricData(g, adjugate(g)/det_g, sqrt(abs(det_g)))
end
function compute_metric_data(md::MetricData{Nothing,<:SymTensorValue,Nothing})
  g_inv = getfield(md, :g_inv)
  det_g = inv(det(g_inv))
  MetricData(adjugate(g_inv)*det_g, g_inv, sqrt(abs(det_g)))
end

"""
    parse_metric_data(metric_args...)
    parse_metric_data(md::MetricData) = md

Gather the provided metric data amongst `g`, `g_inv`, `sqrt_abs_g` given as
positional `metric_args`, in any order, without doing any computation.

```julia
g = SymTensorValue{cov3}(1.:6...)
md = parse_metric_data(g)
md = parse_metric_data(sqrt(abs(det(g))), inv(g)) # same data, g_inv stored instead of g
```

Return a `MetricData` with possibly some `nothing` stored, but computes missing
values if the field is accessed. Since `MetricData` is immutable, the computed
value is returned but not stored, so do not access a field several times.

`LinearAlgebra.I` is parsed into a `MetricData` storing no tensor.
"""
parse_metric_data(md::MetricData) = md
parse_metric_data(I::UniformScaling) = MetricData(I)
function parse_metric_data(metric_args...)
  _parse_metric_args((g=nothing, g_inv=nothing, sqrt_abs_g=nothing), metric_args...)
end
#function parse_metric_data(; g::Union{SymTensorValue{Dg,Tg,L},Nothing}=nothing,
#                             g_inv::Union{SymTensorValue{Di,Ti,L},Nothing}=nothing,
#                             sqrt_abs_g::Ts=nothing, kwargs...) where {Dg,Tg,Di,Ti,Ts,L}
#  MetricData(g, g_inv, sqrt_abs_g)
#end

# The positional arguments are folded one at a time into the named tuple of the
# three data, each assigned by its type. The dispatch and the `D isa Cov` test
# on a type parameter resolve at compile time, so accepting the arguments in any
# order costs nothing at runtime.
_parse_metric_args(nt) = MetricData(nt.g, nt.g_inv, nt.sqrt_abs_g)
_parse_metric_args(nt, a, rest...) = _parse_metric_args(_add_metric_arg(nt, a), rest...)

function _add_metric_arg(nt, a::SymTensorValue{D}) where D
  if D isa Cov
    isnothing(nt.g) || throw(ArgumentError("the metric tensor g was given twice"))
    merge(nt, (g = a,))
  elseif D isa Con
    isnothing(nt.g_inv) || throw(ArgumentError("the inverse metric tensor g_inv was given twice"))
    merge(nt, (g_inv = a,))
  else
    throw(ArgumentError("a metric tensor must carry its variance, Cov for g or Con for g_inv, got $(typeof(a))"))
  end
end
function _add_metric_arg(nt, a::Number)
  isnothing(nt.sqrt_abs_g) || throw(ArgumentError("the scalar sqrt_abs_g was given twice"))
  merge(nt, (sqrt_abs_g = a,))
end
_add_metric_arg(nt, a::MultiValue) = throw(ArgumentError(
  "a metric tensor must be a SymTensorValue with variance, got $(typeof(a))"))
# `nothing` is an absent datum, so that a caller can forward an optional one
# without branching, as in `parse_metric_data(g, sqrt_abs_g_or_nothing)`.
_add_metric_arg(nt, ::Nothing) = nt
_add_metric_arg(nt, a) = throw(ArgumentError(
  "cannot parse metric data from an argument of type $(typeof(a)): expected a MetricData or \
  LinearAlgebra.I on its own, or else the metric tensor g and/or its inverse g_inv as \
  SymTensorValues with variance, optionally with the scalar sqrt_abs_g, in any order"))

function _require_metric_data(a, b)
  throw(ArgumentError(
    "inner(::$(typeof(a)), ::$(typeof(b))) requires metric information, because \
    the tensors carry varianced indices. Pass the metric tensor g and/or its \
    inverse g_inv, or LinearAlgebra.I for the Euclidean inner product, see \
    `parse_metric_data` and `flat_inner`."))
end
_require_metric_data(a, b, metric_args::Vararg{_MetricArg}) = parse_metric_data(metric_args...)

