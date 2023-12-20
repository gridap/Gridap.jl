#############################
# time_derivative interface #
#############################
"""
    time_derivative(f::DerivableType) -> DerivableType

Build the first-order time derivative operator for `f`.
"""
function time_derivative(f)
  @abstractmethod
end

"""
    time_derivative(f::DerivableType, ::Val{k}) -> DerivableType

Build the `k`-th-order time derivative operator for `f`.
"""
function time_derivative(f, ::Val{0})
  f
end

function time_derivative(f, ::Val{1})
  time_derivative(f)
end

function time_derivative(f, ::Val{k}) where {k}
  time_derivative(time_derivative(f), Val(k - 1))
end

"""
    ∂t(f::DerivableType) -> DerivableType

Build the first-th-order time derivative operator for `f`.

Alias for `time_derivative(f)`.
"""
function ∂t(f)
  time_derivative(f)
end

"""
    ∂t(f::DerivableType, ::Val{k}) -> DerivableType

Build the `k`-th-order time derivative operator for `f`.

Alias for `time_derivative(f, Val(k))`.
"""
function ∂t(f, ::Val{k}) where {k}
  time_derivative(f, Val(k))
end

"""
    ∂tt(f::DerivableType) -> DerivableType

Second-order time derivative operator for `f`.

Alias for `time_derivative(f, Val(2))`.
"""
function ∂tt(f)
  time_derivative(f, Val(2))
end

#################################
# Specialisation for `Function` #
#################################
function time_derivative(f::Function)
  function dfdt(x, t)
    z = zero(return_type(f, x, t))
    _time_derivative(f, x, t, z)
  end
  # Extend definition to include restrictions
  _dfdt(x, t) = dfdt(x, t)
  _dfdt(x::VectorValue) = t -> dfdt(x, t)
  _dfdt(t::Real) = x -> dfdt(x, t)
  return _dfdt
end

function _time_derivative(f, x, t, z)
  ForwardDiff.derivative(t -> f(x, t), t)
end

function _time_derivative(f, x, t, z::VectorValue)
  VectorValue(ForwardDiff.derivative(t -> get_array(f(x, t)), t))
  # VectorValue(ForwardDiff.derivative(t -> f(x, t), t))
end

function _time_derivative(f, x, t, z::TensorValue)
  TensorValue(ForwardDiff.derivative(t -> get_array(f(x, t)), t))
end

###############################
# Specialisation for `Number` #
###############################
time_derivative(x::Number) = zero(x)
