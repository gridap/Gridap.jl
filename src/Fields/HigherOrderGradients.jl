#####################################################################################
# Arbitrary-order gradients by forward-mode automatic differentiation
#
# The N-th gradient of a field at a point x is read off the field evaluated at x seeded
# with N nested `ForwardDiff.Dual`s: level L adds a `Dual{NthGradientTag{L}}` with D partials
# (one per coordinate direction) on top of the point of level L-1. The dual type depends on
# (D, N) only, never on the field, so a basis of any length is evaluated ONCE at the seeded
# point, and the derivative compiles once per (basis type, D, N).
#
# This is the default for `FieldGradient{N}` and `FieldGradientArray{N}`: hand-written
# kernels (e.g. N ≤ 2 for `PolynomialBasis`) are more specific and win where they exist.
#
# The cost of a nested dual is (D+1)^N numbers per scalar, so this is meant for moderate N
# (N ≤ 4 in 3D, a bit more in 2D).
#####################################################################################

"""
    struct NthGradientTag{L}

`ForwardDiff` tag of the L-th seeding level of `seed_point`.
"""
struct NthGradientTag{L} end

ForwardDiff.:≺(::Type{NthGradientTag{L1}}, ::Type{NthGradientTag{L2}}) where {L1,L2} = L1 < L2

"""
    seed_point(x::Point{D}, Val(N))

The point `x` seeded with `N` nested duals of width `D`: evaluating a field at it gives a
value whose mixed partials of order `N` are the `N`-th derivatives of the field at `x`.
"""
seed_point(x::Point, ::Val{0}) = x

function seed_point(x::Point{D}, ::Val{N}) where {D,N}
  y = seed_point(x, Val(N-1))
  Point(ntuple(Val(D)) do d
    ForwardDiff.Dual{NthGradientTag{N}}(
      y[d], ForwardDiff.Partials(ntuple(k -> k == d ? one(y[d]) : zero(y[d]), Val(D))))
  end)
end

"""
    nth_gradient_type(::Type{T}, x::Point, Val(N))

Type of the `N`-th gradient of a `T`-valued field, i.e. `gradient_type` applied `N`
times. Scalar values give `VectorValue`, `TensorValue`, `ThirdOrderTensorValue` and then
`HighOrderTensorValue`.
"""
function nth_gradient_type(::Type{T}, x::Point, ::Val{N}) where {T,N}
  G = T
  for _ in 1:N
    G = gradient_type(G, x)
  end
  G
end

_components(v::Number) = (v,)
_components(v::MultiValue) = Tuple(v)
_num_components(::Type{<:Number}) = 1
_num_components(::Type{<:MultiValue{S,T,R,L}}) where {S,T,R,L} = L

"""
    taylor_tensor(::Type{G}, v, Val(N), Val(D))

The tensor of type `G` whose entry `(i₁,…,i_N, c…)` is `∂_{i₁}⋯∂_{i_N} v_c`, read off the value
`v` of a field at an `N`-times seeded point (see `seed_point`). Derivative indices
come first (Gridap's convention, `∇f[i,j] = ∂ᵢf_j`), in column-major order.
"""
@generated function taylor_tensor(::Type{G}, v, ::Val{N}, ::Val{D}) where {G,N,D}
  N == 0 && return :(v)
  ex = Expr[]
  for c in 1:_num_components(v), I in CartesianIndices(ntuple(_ -> D, N))
    e = :(t[$c])
    for k in 1:N
      e = :(ForwardDiff.partials($e, $(I[k])))
    end
    push!(ex, e)
  end
  :(t = _components(v); $G(($(ex...),)))
end

# Single field: the default when no hand-written kernel exists. Not cached: the field is
# evaluated at the seeded point with its own `evaluate`.
function evaluate!(cache, fg::FieldGradient{N}, x::Point{D}) where {N,D}
  fx = evaluate(fg.object, seed_point(x, Val(N)))
  G = nth_gradient_type(typeof(evaluate(fg.object, x)), x, Val(N))
  taylor_tensor(G, fx, Val(N), Val(D))
end

"""
    ad_return_cache(fg::FieldGradientArray{N}, x::AbstractVector{<:Point})
    ad_evaluate!(cache, fg::FieldGradientArray{N}, x::AbstractVector{<:Point})

AD evaluation of the `N`-th gradient of a whole array of fields at a vector of points: the
array is evaluated once at the seeded points. Default of `FieldGradientArray{N}`, and the
fallback of bases that only implement low-order gradients by hand.
"""
function ad_return_cache(fg::FieldGradientArray{N}, x::AbstractVector{<:Point}) where N
  f = fg.fa
  xi = testitem(x)
  xd = CachedArray([seed_point(xj, Val(N)) for xj in x])
  cf = return_cache(f, xd.array)
  G = nth_gradient_type(eltype(return_value(f, x)), xi, Val(N))
  fx = evaluate!(cf, f, xd.array)
  r = CachedArray(zeros(G, size(fx)))
  xd, cf, r
end

function ad_evaluate!(cache, fg::FieldGradientArray{N}, x::AbstractVector{<:Point{D}}) where {N,D}
  xd, cf, r = cache
  setsize!(xd, size(x))
  for i in eachindex(x)
    xd.array[i] = seed_point(x[i], Val(N))
  end
  fx = evaluate!(cf, fg.fa, xd.array)
  setsize!(r, size(fx))
  G = eltype(r.array)
  for i in eachindex(fx)
    @inbounds r.array[i] = taylor_tensor(G, fx[i], Val(N), Val(D))
  end
  r.array
end

function return_cache(fg::FieldGradientArray{N,<:AbstractArray{<:Field}}, x::AbstractVector{<:Point}) where N
  ad_return_cache(fg, x)
end

function evaluate!(cache, fg::FieldGradientArray{N,<:AbstractArray{<:Field}}, x::AbstractVector{<:Point}) where N
  ad_evaluate!(cache, fg, x)
end

# The N-th gradient operator

"""
    struct NthGradient{N} <: Function

The `N`-th gradient operator, `NthGradient{N}()(f) == gradient(f, Val(N))`. Use
`Broadcasting(NthGradient{N}())` on arrays of fields, as for `∇` and `∇∇`.
"""
struct NthGradient{N} <: Function end

(::NthGradient{N})(f) where N = gradient(f, Val(N))

function (g::NthGradient)(a::LinearCombinationField)
  LinearCombinationField(a.values, Broadcasting(g)(a.fields), a.column)
end

evaluate!(cache, g::Broadcasting{<:NthGradient}, a::Field) = g.f(a)

return_value(k::Broadcasting{<:NthGradient}, a::AbstractArray{<:Field}) = evaluate(k, a)

function evaluate!(cache, ::Broadcasting{NthGradient{N}}, a::AbstractArray{<:Field}) where N
  FieldGradientArray{N}(a)
end

function evaluate!(cache, ::Broadcasting{NthGradient{N}}, a::FieldGradientArray{M}) where {N,M}
  FieldGradientArray{N+M}(a.fa)
end

function evaluate!(cache, k::Broadcasting{NthGradient{N}}, a::LinearCombinationFieldVector) where N
  LinearCombinationFieldVector(a.values, k(a.fields))
end

function evaluate!(cache, k::Broadcasting{NthGradient{N}}, a::Transpose{<:Field}) where N
  transpose(k(a.parent))
end

# Same lazy optimisations as for `∇` and `∇∇` (ApplyOptimizations.jl): differentiate the
# reference basis once, not every cell's linear combination.

function lazy_map(k::Broadcasting{<:NthGradient}, a::LazyArray{<:Fill{typeof(linear_combination)}})
  lazy_map(linear_combination, a.args[1], lazy_map(k, a.args[2]))
end

function lazy_map(k::Broadcasting{<:NthGradient}, a::LazyArray{<:Fill{typeof(transpose)}})
  lazy_map(transpose, lazy_map(k, a.args[1]))
end

# Pushforward of the N-th gradient by an affine map
#
# For x = Jξ + b, ∂ₓ = J⁻ᵀ ∂_ξ in every derivative index, so the physical N-th gradient is the
# reference one with each of its first N indices contracted with J⁻ᵀ (`pinvJt`): `push_∇` at
# N = 1 and `push_∇∇` at N = 2.

function _contract_index(a::SArray{S,T}, J, ::Val{k}) where {S,T,k}
  D = size(J, 1)
  b = MArray{S,promote_type(T,eltype(J))}(undef)
  for I in CartesianIndices(a)
    s = zero(eltype(b))
    t = Tuple(I)
    for i in 1:D
      s += J[t[k], i] * a[Base.setindex(t, i, k)...]
    end
    b[I] = s
  end
  SArray(b)
end

"""
    push_∇ⁿ(∇ⁿa::MultiValue, Jt_inv::MultiValue, Val(N))

Pushforward of the reference `N`-th gradient `∇ⁿa` by an affine map with inverse transposed
Jacobian `Jt_inv`: each of the first `N` indices is contracted with `Jt_inv`.
"""
@generated function push_∇ⁿ(t::MultiValue, Jinv::MultiValue, ::Val{N}) where N
  ex = [:(a = _contract_index(a, J, Val($k))) for k in 1:N]
  quote
    a = get_array(t)
    J = get_array(Jinv)
    $(ex...)
    $(Base.typename(t).wrapper){$(t.parameters[1:end-1]...)}(Tuple(a))
  end
end

struct _PushNthGradientMap{N} <: Function end
(::_PushNthGradientMap{N})(t, Jt_inv) where N = push_∇ⁿ(t, Jt_inv, Val(N))

"""
    struct PushNthGradient{N} <: Function

`PushNthGradient{N}()(∇ⁿa, ϕ)` is the pushforward of the reference `N`-th gradient `∇ⁿa` by the
cell map `ϕ`; implemented for affine maps (see `push_∇ⁿ`).
"""
struct PushNthGradient{N} <: Function end

function (::PushNthGradient{N})(∇ⁿa::Field, ϕ::Field) where N
  @notimplemented """\n
  Pushforward of N-th order derivatives of reference quantities is only implemented for
  affine cell maps.
  """
end

function (::PushNthGradient{N})(∇ⁿa::Field, ϕ::AffineField) where N
  Jt_inv = pinvJt(∇(ϕ))
  Operation(_PushNthGradientMap{N}())(∇ⁿa, Jt_inv)
end

function lazy_map(
  k::Broadcasting{PushNthGradient{N}}, cell_∇ⁿa::AbstractArray, cell_map::AbstractArray{<:AffineField}) where N
  cell_invJt = lazy_map(Operation(pinvJt), lazy_map(∇, cell_map))
  lazy_map(Broadcasting(Operation(_PushNthGradientMap{N}())), cell_∇ⁿa, cell_invJt)
end

function lazy_map(
  k::Broadcasting{PushNthGradient{N}},
  cell_∇ⁿat::LazyArray{<:Fill{typeof(transpose)}},
  cell_map::AbstractArray{<:AffineField}) where N
  lazy_map(transpose, lazy_map(k, cell_∇ⁿat.args[1], cell_map))
end
