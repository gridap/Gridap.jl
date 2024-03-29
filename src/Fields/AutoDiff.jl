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
  VectorValue(ForwardDiff.gradient(f,get_array(x)))
end

function gradient(f::Function,x::Point,fx::VectorValue)
  TensorValue(transpose(ForwardDiff.jacobian(y->get_array(f(y)),get_array(x))))
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

function divergence(f::Function,x::Point,fx)
  tr(gradient(f,x,fx))
end

function divergence(f::Function,x::Point,fx::TensorValue{2,2})
  g(x) = SVector(f(x).data)
  a = ForwardDiff.jacobian(g,get_array(x))
  VectorValue(
    a[1,1]+a[2,2],
    a[3,1]+a[4,2],
  )
end

function divergence(f::Function,x::Point,fx::TensorValue{3,3})
  g(x) = SVector(f(x).data)
  a = ForwardDiff.jacobian(g,get_array(x))
  VectorValue(
    a[1,1]+a[2,2]+a[3,3],
    a[4,1]+a[5,2]+a[6,3],
    a[7,1]+a[8,2]+a[9,3],
   )
end

function divergence(f::Function,x::Point,fx::SymTensorValue{2})
  g(x) = SVector(f(x).data)
  a = ForwardDiff.jacobian(g,get_array(x))
  VectorValue(
    a[1,1]+a[2,2],
    a[2,1]+a[3,2],
  )
end

function divergence(f::Function,x::Point,fx::SymTensorValue{3})
  g(x) = SVector(f(x).data)
  a = ForwardDiff.jacobian(g,get_array(x))
  VectorValue(
    a[1,1]+a[2,2]+a[3,3],
    a[2,1]+a[4,2]+a[5,3],
    a[3,1]+a[5,2]+a[6,3],
   )
end

function curl(f::Function,x::Point)
  grad2curl(gradient(f,x))
end

function laplacian(f::Function,x::Point)
  laplacian(f,x,return_value(f,x))
end

function laplacian(f::Function,x::Point,fx::Number)
  tr(ForwardDiff.jacobian(y->ForwardDiff.gradient(f,y), get_array(x)))
end

function laplacian(f::Function,x::Point,fx::VectorValue)
  A = length(x)
  B = length(fx)
  a = ForwardDiff.jacobian(y->transpose(ForwardDiff.jacobian(z->get_array(f(z)),y)), get_array(x))
  tr(ThirdOrderTensorValue{A,A,B}(Tuple(transpose(a))))
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
