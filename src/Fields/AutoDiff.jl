# Number differentiation

function gradient(f::Number)
  function grad_f(x::Point)
    zero(return_type(outer,x,f))
  end
end

function ∇∇(f::Number)
  function hess_f(x::Point)
    g = gradient(f)(x)
    gradient(g)(x)
  end
end

# Automatic differentiation of functions

function gradient(f::Function)
  function _gradient(x)
    gradient(f,x)
  end
end

function symmetric_gradient(f::Function)
  function _symmetric_gradient(x)
    symmetric_gradient(f,x)
  end
end

function divergence(f::Function)
  function _divergence(x)
    divergence(f,x)
  end
end

function curl(f::Function)
  function _curl(x)
    curl(f,x)
  end
end

function laplacian(f::Function)
  function _laplacian(x)
    laplacian(f,x)
  end
end

# Specialization when x::Point

function gradient(f::Function,x::Point)
  gradient(f,x,return_value(f,x))
end

function gradient(f::Function,x::Point,fx::Number)
  VectorValue(ForwardDiff.gradient(f∘Point,get_array(x)))
end

# Implementation for tensor values of any order
# Does not exploit possible symmetries
function gradient(f::Function,x::Point,fx::MultiValue)
  # Jacobian of size (num_components(fx), D), components in column-major order. Transposing
  # puts the derivative index first, which is the storage order of `gradient_type`.
  a = ForwardDiff.jacobian(get_array∘f∘Point,get_array(x))
  G = gradient_type(typeof(fx),x)
  G(Tuple(transpose(a)))
end

function symmetric_gradient(f::Function,x::Point)
  symmetric_part(gradient(f,x))
end

function divergence(f::Function,x::Point)
  divergence(f,x,return_value(f,x))
end

function divergence(f::Function,x::Point,fx::VectorValue)
  tr(gradient(f,x,fx))
end

function divergence(f::Function,x::Point{D},fx::S) where S<:MultiValue{Tuple{D,A},T} where {D,A,T}
  a = ForwardDiff.jacobian(get_array∘f∘Point,get_array(x))
  VectorValue{A,T}( ntuple(k -> sum(i-> a[(k-1)*D+i,i], 1:D),A) )
end

function divergence(f::Function,x::Point{D},fx::S) where S<:MultiValue{Tuple{D,A,B},T} where {D,A,B,T}
  a = ForwardDiff.jacobian(get_array∘f∘Point,get_array(x))
  TensorValue{A,B,T}( ntuple(k -> sum(i-> a[(k-1)*D+i,i], 1:D),A*B) )
end

function divergence(f::Function,x::Point,fx::TensorValue{2,2})
  g(x) = SVector(f(Point(x)).data)
  a = ForwardDiff.jacobian(g,get_array(x))
  VectorValue(
    a[1,1]+a[2,2],
    a[3,1]+a[4,2],
  )
end

function divergence(f::Function,x::Point,fx::TensorValue{3,3})
  g(x) = SVector(f(Point(x)).data)
  a = ForwardDiff.jacobian(g,get_array(x))
  VectorValue(
    a[1,1]+a[2,2]+a[3,3],
    a[4,1]+a[5,2]+a[6,3],
    a[7,1]+a[8,2]+a[9,3],
   )
end

function divergence(f::Function,x::Point{2},fx::AbstractSymTensorValue{2})
  g(x) = SVector(f(Point(x)).data)
  a = ForwardDiff.jacobian(g,get_array(x))
  VectorValue(
    a[1,1]+a[2,2],
    a[2,1]+a[3,2],
  )
end

function divergence(f::Function,x::Point{3},fx::AbstractSymTensorValue{3})
  g(x) = SVector(f(Point(x)).data)
  a = ForwardDiff.jacobian(g,get_array(x))
  VectorValue(
    a[1,1]+a[2,2]+a[3,3],
    a[2,1]+a[4,2]+a[5,3],
    a[3,1]+a[5,2]+a[6,3],
   )
end

function divergence(f::Function,x::Point,fx::MultiValue)
  @notimplemented
end

function curl(f::Function,x::Point)
  grad2curl(gradient(f,x))
end

function laplacian(f::Function,x::Point)
  laplacian(f,x,return_value(f,x))
end

function laplacian(f::Function,x::Point,fx::Number)
  tr(ForwardDiff.jacobian(y->ForwardDiff.gradient(f∘Point,y), get_array(x)))
end

function laplacian(f::Function,x::Point{A},fx::VectorValue{B,T}) where {A,B,T}
  a = MMatrix{A*B,A,T}(undef)
  ForwardDiff.jacobian!(a, y->ForwardDiff.jacobian(get_array∘f∘Point,y), get_array(x))
  VectorValue{B,T}( ntuple(k -> sum(i-> a[(i-1)*B+k,i], 1:A),B) )
end

function laplacian(f::Function,x::Point{A},fx::S) where S<:MultiValue{Tuple{B,C},T} where {A,B,C,T}
  a = MMatrix{A*B*C,A,T}(undef)
  ForwardDiff.jacobian!(a, y->ForwardDiff.jacobian(get_array∘f∘Point,y), get_array(x))
  t = ntuple(k -> sum(i-> a[(i-1)*B*C+k,i], 1:A),B*C)
  S(SMatrix{B,C}(t)) #Necessary cast to build e.g. symmetric tensor values
end

# Implementation for any third order tensor values
function laplacian(f::Function,x::Point{A},fx::S) where S<:MultiValue{Tuple{B,C,D},T} where {A,B,C,D,T}
  a = MMatrix{A*B*C*D,A,T}(undef)
  ForwardDiff.jacobian!(a, y->ForwardDiff.jacobian(get_array∘f∘Point,y), get_array(x))
  t = ntuple(k -> sum(i-> a[(i-1)*B*C*D+k,i], 1:A),B*C*D)
  S(SArray{Tuple{B,C,D}}(t))
end

function laplacian(f::Function,x::Point,fx::MultiValue)
  @notimplemented
end

# Specialization when x::StaticVector

gradient(f::Function,x::SVector) = gradient(f,Point(x))
symmetric_gradient(f::Function,x::SVector) = symmetric_gradient(f,Point(x))
divergence(f::Function,x::SVector) = divergence(f,Point(x))
curl(f::Function,x::SVector) = curl(f,Point(x))
laplacian(f::Function,x::SVector) = laplacian(f,Point(x))


##############################################
# Variance consistent differential operators #
##############################################

# `N`th order partial derivatives of an `L` component function `g`. Each
# nesting flattens the previous result column major, so with `D = length(y)`:
#     a[l + L*∑_{k<N} (i_k-1)*D^(k-1), i_N] == ∂_{i_1} ⋯ ∂_{i_N} g_l
_nested_jacobian(g, y::SVector{D,T}, ::Val{1}, ::Val{L}) where {D,T,L} =
  ForwardDiff.jacobian(g,y)::SMatrix{L,D,T,L*D}

function _nested_jacobian(g, y::SVector{D,T}, ::Val{N}, ::Val{L}) where {D,T,N,L}
  ForwardDiff.jacobian(z -> _nested_jacobian(g,z,Val(N-1),Val(L)), y)::SMatrix{L*D^(N-1),D,T,L*D^N}
end

@generated function _partials_data(a::StaticMatrix, ::Val{D}, ::Val{N}, ::Val{L}) where {D,N,L}
  # A row of `a` holds one component and the `N-1` outermost derivation indices,
  # the innermost one indexing the columns.
  rows = LinearIndices((L, ntuple(_->D, N-1)...))
  comps = [ :( a[$(rows[l, I[1:N-1]...]), $(I[N])] )
            for l in 1:L for I in sorted_multisets(D,N) ]
  Expr(:tuple, comps...)
end

function (k::PartialDerivatives)(f::Function)
  function _partial_derivatives(x)
    k(f,x)
  end
end

function (k::PartialDerivatives)(f::Function, x::Point)
  k(f, x, return_value(f,x))
end

(::PartialDerivatives{0})(f::Function, x::Point, fx::Number) = fx*one(eltype(x))

# Only the independent components of `fx` are differentiated, and only the
# non-decreasing derivation multi-indices are stored in the resulting ProductTensorValue.
function (::PartialDerivatives{N})(f::Function, x::Point{D}, fx::Number) where {D,N}
  V = derivatives_type(typeof(fx), x, Val(N))
  L = num_indep_components(typeof(fx))
  iszero(L) && return zero(V)
  a = _nested_jacobian(_indep_comp_svector∘f∘Point, get_array(x), Val(N), Val(L))
  V(_partials_data(a, Val(D), Val(N), Val(L)))
end

_indep_comp_svector(a::Number) = SVector(get_indep_components(a))

function exterior_derivative(ω::Function)
  function _exterior_derivative(x)
    exterior_derivative(ω,x)
  end
end

exterior_derivative(ω::Function, x::Point) =
  partials_to_exteriorderivative(partial_derivatives(ω,x))

function codifferential(ω::Function, I::UniformScaling)
  function _codifferential(x)
    codifferential(ω,x,I)
  end
end

codifferential(ω::Function, x::Point, ::UniformScaling) =
  partials_to_flat_codifferential(partial_derivatives(ω,x))

# With the metric given as functions too, `g` and the `metric_args` that
# `partials_to_codifferential` takes after it:
#
# function codifferential(ω::Function, g::Function, metric_args::Function...)
#   function _codifferential(x)
#     codifferential(ω,x,g,metric_args...)
#   end
# end
#
# function codifferential(ω::Function, x::Point, g::Function, metric_args::Function...)
#   gx = g(x)
#   args = map(a -> a(x), metric_args)
#   Γ = partials_to_Christoffel_2(partial_derivatives(g,x), gx, args...)
#   partials_to_codifferential(partial_derivatives(ω,x), ω(x), Γ, gx, args...)
# end

function lie_derivative(v::Function, ω::Function)
  function _lie_derivative(x)
    lie_derivative(v,ω,x)
  end
end

lie_derivative(v::Function, ω::Function, x::Point) = partials_to_liederivative(
  v(x), ω(x), partial_derivatives(v,x), partial_derivatives(ω,x))

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

struct NthGradientTag{L} end

ForwardDiff.:≺(::Type{NthGradientTag{L1}}, ::Type{NthGradientTag{L2}}) where {L1,L2} = L1 < L2

seed_point(x::Point, ::Val{0}) = x

function seed_point(x::Point{D}, ::Val{N}) where {D,N}
  y = seed_point(x, Val(N-1))
  Point(ntuple(Val(D)) do d
    ForwardDiff.Dual{NthGradientTag{N}}(
      y[d], ForwardDiff.Partials(ntuple(k -> k == d ? one(y[d]) : zero(y[d]), Val(D))))
  end)
end

#     taylor_tensor(::Type{G}, v, Val(N), Val(D))
#
# The tensor of type `G` whose entry `(i₁,…,i_N, c…)` is `∂_{i₁}⋯∂_{i_N} v_c`, read off the value
# `v` of a field at an `N`-times seeded point. Derivative indices
# come first (Gridap's convention, `∇f[i,j] = ∂ᵢf_j`), in column-major order.
@generated function taylor_tensor(::Type{G}, v, ::Val{N}, ::Val{D}) where {G,N,D}
  N == 0 && return :(v)
  ex = Expr[]
  for c in 1:num_components(v), I in CartesianIndices(ntuple(_ -> D, N))
    e = :(v[$c])
    for k in 1:N
      e = :(ForwardDiff.partials($e, $(I[k])))
    end
    push!(ex, e)
  end
  return :($G(($(ex...),)))
end
