
"""
    divergence(f)
"""
divergence(f) = tr(gradient(f))

function symmetric_gradient end

"""
    symmetric_gradient(f)
"""
symmetric_gradient(f) = symmetric_part(gradient(f))

"""
    const ε = symmetric_gradient

Alias for the symmetric gradient
"""
const ε = symmetric_gradient

"""
    curl(f)
"""
curl(f) = grad2curl(gradient(f))

"""
    grad2curl(∇f)
"""
function grad2curl(f)
  @abstractmethod
end

function laplacian end

"""
    const Δ = laplacian

Alias for the `laplacian` function
"""
const Δ = laplacian

"""
    laplacian(f)
"""
function laplacian(f)
  g = gradient(f)
  h = gradient(g)
  tr(h)
end

"""
    ∇⋅f

Equivalent to

    divergence(f)
"""
dot(::typeof(∇),f) = divergence(f)
dot(::typeof(∇),f::GridapType) = divergence(f)

function (*)(::typeof(∇),f)
  msg = "Syntax ∇*f has been removed, use ∇⋅f (\\nabla \\cdot f) instead"
  error(msg)
end

function (*)(::typeof(∇),f::GridapType)
  msg = "Syntax ∇*f has been removed, use ∇⋅f (\\nabla \\cdot f) instead"
  error(msg)
end

"""
    outer(∇,f)

Equivalent to

    gradient(f)
"""
outer(::typeof(∇),f) = gradient(f)
outer(::typeof(∇),f::GridapType) = gradient(f)

"""
    outer(f,∇)

Equivalent to

    transpose(gradient(f))
"""
outer(f,::typeof(∇)) = transpose(gradient(f))
outer(f::GridapType,::typeof(∇)) = transpose(gradient(f))

"""
    cross(∇,f)

Equivalent to
    
    curl(f)
"""
cross(::typeof(∇),f) = curl(f)
cross(::typeof(∇),f::GridapType) = curl(f)

# Helpers

grad2curl(f::Field) = apply_kernel_to_field(bcast(_curl_kernel),f)

grad2curl(f::AbstractArray{<:Field}) = apply_to_field_array(bcast(_curl_kernel),f)

grad2curl(::Type{T}, f::AbstractArray{<:Field}) where T = apply_to_field_array(T,bcast(_curl_kernel),f)

grad2curl(∇u::MultiValue) = _curl_kernel(∇u)

function _curl_kernel(∇u::TensorValue{2})
  ∇u[1,2] - ∇u[2,1]
end

function _curl_kernel(∇u::TensorValue{3})
  c1 = ∇u[2,3] - ∇u[3,2]
  c2 = ∇u[3,1] - ∇u[1,3]
  c3 = ∇u[1,2] - ∇u[2,1]
  VectorValue(c1,c2,c3)
end

# Automatic differentiation of functions

function gradient(f::Function)
  function grad_f(x)
    _grad_f(f,x,zero(return_type(f,typeof(x))))
  end
end

function _grad_f(f,x,fx)
  VectorValue(ForwardDiff.gradient(f,get_array(x)))
end

function _grad_f(f,x,fx::VectorValue)
  TensorValue(transpose(ForwardDiff.jacobian(y->get_array(f(y)),get_array(x))))
end

function _grad_f(f,x,fx::MultiValue)
  @notimplemented
end

function divergence(f::Function)
  x -> tr(ForwardDiff.jacobian(y->get_array(f(y)),get_array(x)))
end

function curl(f::Function)
  x -> grad2curl(TensorValue(transpose(ForwardDiff.jacobian(y->get_array(f(y)),get_array(x)))))
end

function laplacian(f::Function)
  function lapl_f(x)
    _lapl_f(f,x,zero(return_type(f,typeof(x))))
  end
end

function _lapl_f(f,x,fx)
  tr(ForwardDiff.jacobian(y->ForwardDiff.gradient(f,y), get_array(x)))
end

function _lapl_f(f,x,fx::VectorValue)
  A = length(x)
  B = length(fx)
  a = ForwardDiff.jacobian(y->transpose(ForwardDiff.jacobian(z->get_array(f(z)),y)), get_array(x))
  tr(ThirdOrderTensorValue{A,A,B}(Tuple(transpose(a))))
end

function _lapl_f(f,x,fx::MultiValue)
  @notimplemented
end

function symmetric_gradient(f::Function)
    x -> symmetric_part(_grad_f(f,x,zero(return_type(f,typeof(x)))))
end

