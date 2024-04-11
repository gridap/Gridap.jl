#####################
# TimeSpaceFunction #
#####################
"""
    struct TimeSpaceFunction{F} <: Function end

`TimeSpaceFunction` allows for convenient ways to apply differential operators to
functions that depend on time and space. More precisely, if `f` is a function that, to
a given time, returns a function of space (i.e. `f` is evaluated at time `t` and position
`x` via `f(t)(x)`), then `F = TimeSpaceFunction(f)` supports the following syntax:
* `op(F)`: a `TimeSpaceFunction` representing both `t -> x -> op(f)(t)(x)` and `(t, x) -> op(f)(t)(x)`,
* `op(F)(t)`: a function of space representing `x -> op(f)(t)(x)`
* `op(F)(t, x)`: the quantity `op(f)(t)(x)` (this notation is equivalent to `op(F)(t)(x)`),
for all spatial and temporal differential operator, i.e. `op` in `(time_derivative,
gradient, symmetric_gradient, divergence, curl, laplacian)` and their symbolic aliases
(`∂t`, `∂tt`, `∇`, ...).
"""
struct TimeSpaceFunction{F} <: Function
  f::F
end

(ts::TimeSpaceFunction)(t) = ts.f(t)
(ts::TimeSpaceFunction)(t, x) = ts.f(t)(x)

##################################
# Spatial differential operators #
##################################
# Using the rule spatial_op(ts)(t)(x) = spatial_op(ts(t))(x)
for op in (:(Fields.gradient), :(Fields.symmetric_gradient), :(Fields.divergence),
  :(Fields.curl), :(Fields.laplacian))
  @eval begin
    function ($op)(ts::TimeSpaceFunction)
      function _op(t)
        _op_t(x) = $op(ts(t))(x)
        _op_t
      end
      TimeSpaceFunction(_op)
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
  function dfdt(t)
    ft = f(t)
    function dfdt_t(x)
      T = return_type(ft, x)
      _time_derivative(T, f, t, x)
    end
    dfdt_t
  end
  dfdt
end

function _time_derivative(T::Type{<:Real}, f, t, x)
  partial(t) = f(t)(x)
  ForwardDiff.derivative(partial, t)
end

function _time_derivative(T::Type{<:VectorValue}, f, t, x)
  partial(t) = get_array(f(t)(x))
  VectorValue(ForwardDiff.derivative(partial, t))
end

function _time_derivative(T::Type{<:TensorValue}, f, t, x)
  partial(t) = get_array(f(t)(x))
  TensorValue(ForwardDiff.derivative(partial, t))
end

##########################################
# Specialisation for `TimeSpaceFunction` #
##########################################
function time_derivative(ts::TimeSpaceFunction)
  TimeSpaceFunction(time_derivative(ts.f))
end

###############################
# Specialisation for `Number` #
###############################
time_derivative(x::Number) = zero(x)
