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

function gradient(f::Function,x::Point,fx::VectorValue)
  TensorValue(transpose(ForwardDiff.jacobian(get_array∘f∘Point,get_array(x))))
end

# Implementation for all second order tensor values
# Does not exploit possible symmetries
function gradient(f::Function,x::Point{A},fx::S) where S<:MultiValue{Tuple{B,C}} where {A,B,C}
  a = transpose(ForwardDiff.jacobian(get_array∘f∘Point,get_array(x)))
  ThirdOrderTensorValue(reshape(a, (A,B,C)))
end

function gradient(f::Function,x::Point,fx::MultiValue)
  @notimplemented
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

function partial_derivatives(f::Function, ::Val{N}) where N
  function _partial_derivatives(x)
    partial_derivatives(f,x,Val(N))
  end
end

function partial_derivatives(f::Function, x::Point, ::Val{N}) where N
  partial_derivatives(f, x, return_value(f,x), Val(N))
end

partial_derivatives(f::Function, x::Point, fx::Number, ::Val{0}) = fx*one(eltype(x))

# Only the independent components of `fx` are differentiated, and only the
# non-decreasing derivation multi-indices are stored in the resulting ProductTensorValue.
function partial_derivatives(f::Function, x::Point{D}, fx::Number, ::Val{N}) where {D,N}
  V = derivatives_type(typeof(fx), x, Val(N))
  L = num_indep_components(typeof(fx))
  iszero(L) && return zero(V)
  a = _nested_jacobian(_indep_comp_svector∘f∘Point, get_array(x), Val(N), Val(L))
  V(_partials_data(a, Val(D), Val(N), Val(L)))
end

_indep_comp_svector(a::Number) = SVector(get_indep_components(a))
