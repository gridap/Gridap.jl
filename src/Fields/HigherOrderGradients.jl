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

"""
    struct PushNthGradient{N} <: Function

`PushNthGradient{N}()(∇ⁿa, ϕ)` is the pushforward of the reference `N`-th gradient `∇ⁿa` by the
cell map `ϕ`; implemented for affine maps (see `push_∇ⁿ`).
"""
struct PushNthGradient{N} <: Function end

(::PushNthGradient{N})(∇ⁿa::Field, ϕ::Field) where N = push_∇ⁿ(∇ⁿa, ϕ, Val(N))

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
