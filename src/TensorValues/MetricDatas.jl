# Metric data handler

"""
    struct MetricData{
        G    <: Union{SymTensorValue,Nothing},
        Ginv <: Union{SymTensorValue,Nothing},
        T    <: Union{Number,Nothing}
    }
      g::G
      g_inv::Ginv
      sqrt_abs_g::T
    end

Helper struct to handle the metric information consisting of the metric tensor
`g`, its inverse `g_inv`, and `sqrt_abs_g`=√|det(g)|. At least one of the
tensors needs be provided. Missing field are computed on the fly during access.

See [`parse_metric_data`](@ref) for lazy computation, and
[`compute_metric_data`](@ref) for greedy computation.

```
_g = SymTensorValue{cov2}(1.,2,3)
_g_inv = inv(g)
_sqrt_abs_g = sqrt(abs(det(_g)))

g_or_ginv      = rand(Bool) ? _g : _g_inv
det_or_nothing = rand(Bool) ? _sqrt_abs_g : nothing
md = compute_metric_data(g_or_ginv, det_or_nothing)
md.g ≈ _g && md.g_inv ≈ _g_inv && md.sqrt_abs_g ≈ _sqrt_abs_g # true
(; g, g_inv, sqrt_abs_g) = md
   g ≈ _g &&    g_inv ≈ _g_inv &&    sqrt_abs_g ≈ _sqrt_abs_g  # true
```
"""
struct MetricData{
  G    <: Union{SymTensorValue,Nothing},
  Ginv <: Union{SymTensorValue,Nothing},
  T    <: Union{Number,Nothing}
}
  g::G
  g_inv::Ginv
  sqrt_abs_g::T

  function MetricData(g::Union{SymTensorValue{Dcov,Tg,L},Nothing},
                      g_inv::Union{SymTensorValue{Dcon,Ti,L},Nothing},
                      sq_g::Ts) where {Dcov,Tg,Dcon,Ti,L,Ts}

    isnothing(g) && isnothing(g_inv) && throw(ArgumentError("either the g or g_inv kwarg must be provided, got neither"))

    # validate dimension and variance
    _Dcov = isnothing(g)     ? inverse_variance(Dcon) : Dcov
    _Dcon = isnothing(g_inv) ? inverse_variance(Dcov) : Dcon
    @assert _Dcov == _Dcon "metric tensor g and inverse g_inv do not have the same dimension"
    @assert _Dcov isa Cov  "metric tensor g must be covariant"
    @assert _Dcon isa Con  "metric tensor g must be covariant"

    # promote to same scalar type
    _Tg = isnothing(g)     ? Bool : Tg
    _Ti = isnothing(g_inv) ? Bool : Ti
    _Ts = isnothing(sq_g)  ? Bool : Ts
    _T = Base.promote_op(inv, promote_type(_Tg, _Ti, _Ts))
    _g     = isnothing(g)     ? g     : convert(SymTensorValue{_Dcov,_T}, g)
    _g_inv = isnothing(g_inv) ? g_inv : convert(SymTensorValue{_Dcon,_T}, g_inv)
    _sq_g  = isnothing(sq_g)  ? sq_g  : convert(_T, sq_g)

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
    !isnothing(g) && return sqrt(abs(det(g)))
    g_inv = getfield(md, :g_inv)
    return inv(sqrt(abs(det(g))))

  else
    throw(Base.FieldError(MetricData, name))
  end
end

#    compute_metric_data(; g=nothing, g_inv=nothing, sqrt_abs_g=nothing)
#
#
#The `MetricData` constructors will store all three data, eagrly computing all
#those missing (`nothing`) in the argument. This is meant for use at the call site:
#
#```julia
#md = MetricData(g) # or g_inv
#md = MetricData(; g_inv=inv(g), sqrt_abs_g=sqrt(abs(det(g))) )
#```
#
#This is useful if the caller does not now if only either `g` and `g_inv` might
#be needed. Otherwise, the caller may only pass one of them, and let the
#operator compute missing data.
#
"""
    compute_metric_data(g_or_ginv::SymTensorValue)

The callee will use [`parse_metric_data(; kwargs...)`](@ref parse_metric_data)
which lazily compute missing required data.
"""
function compute_metric_data(g_or_ginv::SymTensorValue{D}, sqrt_abs_g=nothing) where {D}
  data = if D isa Con
    # g_or_ginv is g_inv
    det_g = inv(det(g_or_ginv))
    adj_ginv = adjugate(g_or_ginv)
    inverse_variance(D), D, adj_ginv*det_g, g_or_ginv, det_g
  elseif D isa Cov
    # g_or_ginv is g
    det_g = det(g_or_ginv)
    g_adj = adjugate(g_or_ginv)
    D, inverse_variance(D), g_or_ginv, g_adj/det_g, det_g
  else
    throw(ArgumentError("the argument has no variance, got size $(size(g_or_ginv))"))
  end
  Dcov, Dcon, g, g_inv, det_g = data
  sqrt_abs_g = sqrt(abs(det_g))
  MetricData(g, g_inv, sqrt_abs_g)
end

"""
    parse_metric_data(g_or_ginv::SymTensorValue, sqrt_abs_g=nothing)

Extract the available metric data from the given keyword argument, without
doing any computation. Return a `MetricData` with possibly some `nothing`
stored, but computes missing values if the field is accessed. Since
`MetricData` is immutable, the computed value is returned but not stored, so do
not access a field several times.

```julia
g = SymTensorValue{cov3}(1.:6...)
md = parse_metric_data(; g)
```
"""
function parse_metric_data(g_or_ginv::SymTensorValue{Dg}, sqrt_abs_g=nothing) where Dg
  Dg isa Cov && return MetricData(g_or_ginv, nothing, sqrt_abs_g)
  Dg isa Con && return MetricData(nothing, g_or_ginv, sqrt_abs_g)
  throw(ArgumentError("The given (inverse) metric tensor has no variance"))
end
function parse_metric_data(; g::Union{SymTensorValue{Dg,Tg,L},Nothing}=nothing,
                             g_inv::Union{SymTensorValue{Di,Ti,L},Nothing}=nothing,
                             sqrt_abs_g::Ts=nothing, kwargs...) where {Dg,Tg,Di,Ti,Ts,L}
  MetricData(g, g_inv, sqrt_abs_g)
end


