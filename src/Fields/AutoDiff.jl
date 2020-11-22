# Number differentiation

function gradient(f::Number)
  @inline function grad_f(x::Point)
    zero(return_type(outer,x,f))
  end
end

function ∇∇(f::Number)
  @inline function hess_f(x::Point)
    g = gradient(f)(x)
    gradient(g)(x)
  end
end

# Automatic differentiation of functions

function gradient(f::Function)
  @inline function _gradient(x)
    gradient(f,x)
  end
end

function symmetric_gradient(f::Function)
  @inline function _symmetric_gradient(x)
    symmetric_gradient(f,x)
  end
end

function divergence(f::Function)
  @inline function _divergence(x)
    divergence(f,x)
  end
end

function curl(f::Function)
  @inline function _curl(x)
    curl(f,x)
  end
end

function laplacian(f::Function)
  @inline function _laplacian(x)
    laplacian(f,x)
  end
end

# Specialization when x::Point

function gradient(f::Function,x::Point)
  gradient(f,x,return_value(f,x))
end

@inline function gradient(f::Function,x::Point,fx::Number)
  VectorValue(ForwardDiff.gradient(f,get_array(x)))
end

@inline function gradient(f::Function,x::Point,fx::VectorValue)
  TensorValue(transpose(ForwardDiff.jacobian(y->get_array(f(y)),get_array(x))))
end

function gradient(f::Function,x::Point,fx::MultiValue)
  @notimplemented
end

function symmetric_gradient(f::Function,x::Point)
  symmetric_part(gradient(f,x))
end

function divergence(f::Function,x::Point)
  tr(gradient(f,x))
end

function curl(f::Function,x::Point)
  grad2curl(gradient(f,x))
end

function laplacian(f::Function,x::Point)
  laplacian(f,x,return_value(f,x))
end

@inline function laplacian(f::Function,x::Point,fx::Number)
  tr(ForwardDiff.jacobian(y->ForwardDiff.gradient(f,y), get_array(x)))
end

@inline function laplacian(f::Function,x::Point,fx::VectorValue)
  A = length(x)
  B = length(fx)
  a = ForwardDiff.jacobian(y->transpose(ForwardDiff.jacobian(z->get_array(f(z)),y)), get_array(x))
  tr(ThirdOrderTensorValue{A,A,B}(Tuple(transpose(a))))
end

function laplacian(f::Function,x::Point,fx::MultiValue)
  @notimplemented
end

# Specialization when x::StaticVector

@inline gradient(f::Function,x::SVector) = gradient(f,Point(x))
@inline symmetric_gradient(f::Function,x::SVector) = symmetric_gradient(f,Point(x))
@inline divergence(f::Function,x::SVector) = divergence(f,Point(x))
@inline curl(f::Function,x::SVector) = curl(f,Point(x))
@inline laplacian(f::Function,x::SVector) = laplacian(f,Point(x))
