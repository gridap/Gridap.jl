################
# Time slicing #
################
"""
    struct TimeSlicing{F} <: Function end

`TimeSlicing` is an operator that can be applied to a function of time and space, and
that, to a given time, produces the function of space at that time, allowing for the
following syntax: if `f(t, x)` is defined, then `TimeSlicing(f)(t)(x) = f(t, x)`.

This is needed prevent performance drops with automatic differentiation (e.g. when taking
the time derivative of Dirichlet boundary conditions), and with integration (i.e. when
an integrand contains a time-dependent coefficient).
"""
struct TimeSlicing{F} <: Function
  f::F
end

function (ts::TimeSlicing)(t)
  _ftx(x) = ts.f(t, x)
  _ftx
end

function (ts::TimeSlicing)(t, x)
  ts.f(t, x)
end

"""
    const time_slicing = TimeSlicing

Alias for `TimeSlicing`.
"""
const time_slicing = TimeSlicing

##################################
# Spatial differential operators #
##################################
# Using the rule spatial_op(ts)(t, x) = spatial_op(ts(t))(x)
for op in (:(Fields.gradient), :(Fields.symmetric_gradient), :(Fields.divergence),
  :(Fields.curl), :(Fields.laplacian))
  @eval begin
    function ($op)(ts::TimeSlicing)
      _op(t, x) = $op(ts(t))(x)
      TimeSlicing(_op)
    end
  end
end

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
  function dfdt(t, x)
    T = return_type(f, t, x)
    _time_derivative(T, f, t, x)
  end
  dfdt
end

function _time_derivative(T::Type{<:Real}, f, t, x)
  partial(t) = f(t, x)
  ForwardDiff.derivative(partial, t)
end

function _time_derivative(T::Type{<:VectorValue}, f, t, x)
  partial(t) = get_array(f(t, x))
  VectorValue(ForwardDiff.derivative(partial, t))
end

function _time_derivative(T::Type{<:TensorValue}, f, t, x)
  partial(t) = get_array(f(t, x))
  TensorValue(ForwardDiff.derivative(partial, t))
end

####################################
# Specialisation for `TimeSlicing` #
####################################
function time_derivative(ts::TimeSlicing)
  TimeSlicing(time_derivative(ts.f))
end

###############################
# Specialisation for `Number` #
###############################
time_derivative(x::Number) = zero(x)
