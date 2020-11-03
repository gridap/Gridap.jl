
"""
    abstract type FEOperator <: GridapType

A `FEOperator` contains finite element problem,
that is assembled as far as possible and ready to be solved.
See also [FETerm](@ref)
"""
abstract type FEOperator <: GridapType end

"""
"""
function get_test(op::FEOperator)
  @abstractmethod
end

"""
"""
function get_trial(op::FEOperator)
  @abstractmethod
end

"""
"""
function allocate_residual(op::FEOperator,u::FEFunction)
  @abstractmethod
end

"""
    $(SIGNATURES)

Inplace version of [`residual`](@ref).
"""
function residual!(b::AbstractVector,op::FEOperator,u::FEFunction)
  @abstractmethod
end

"""
    $(SIGNATURES)

Compute the residual of `op` at `u`. See also [`residual_and_jacobian`](@ref)
"""
function residual(op::FEOperator,u::FEFunction)
  b = allocate_residual(op,u)
  residual!(b,op,u)
  b
end

"""
"""
function allocate_jacobian(op::FEOperator,u::FEFunction)
  @abstractmethod
end

"""
    $(SIGNATURES)

Inplace version of [`jacobian`](@ref).
"""
function jacobian!(A::AbstractMatrix,op::FEOperator,u::FEFunction)
  @abstractmethod
end

"""
    $(SIGNATURES)

Compute the jacobian of an operator `op`.
See also [`get_algebraic_operator`](@ref), [`residual_and_jacobian!`](@ref).
"""
function jacobian(op::FEOperator,u::FEFunction)
  A = allocate_jacobian(op,u)
  jacobian!(A,op,u)
  A
end

"""
    $(SIGNATURES)

Inplace version of [`residual_and_jacobian`](@ref).
"""
function residual_and_jacobian!(b::AbstractVector,A::AbstractMatrix,op::FEOperator,u::FEFunction)
  residual!(b,op,u)
  jacobian!(A,op,u)
  (b,A)
end

"""
    residual, jacobian = $(SIGNATURES)

Compute the residual and jacobian of an operator `op` at a given point `u`.
Depending on the nature of `op` the point `u` can either be a plain array or a `FEFunction`.

See also [`jacobian`](@ref), [`residual`](@ref), [`get_algebraic_operator`](@ref).
"""
function residual_and_jacobian(op::FEOperator,u::FEFunction)
  b = residual(op,u)
  A = jacobian(op,u)
  (b,A)
end

"""
"""
function test_fe_operator(op::FEOperator,args...;kwargs...)
  test = get_test(op)
  trial = get_trial(op)
  @test isa(test,FESpace)
  @test isa(trial,FESpace)
  u = zero(trial)
  b = allocate_residual(op,u)
  @test isa(b,AbstractVector)
  residual!(b,op,u)
  b2 = residual(op,u)
  @test isa(b2,AbstractVector)
  A = allocate_jacobian(op,u)
  @test isa(A,AbstractMatrix)
  jacobian!(A,op,u)
  A2 = jacobian(op,u)
  @test isa(A2,AbstractMatrix)
  residual_and_jacobian!(b,A,op,u)
  b, A = residual_and_jacobian(op,u)
  @test isa(b,AbstractVector)
  @test isa(A,AbstractMatrix)
  _op = get_algebraic_operator(op)
  test_nonlinear_operator(_op,args...;kwargs...)
end

"""
    $(SIGNATURES)

Return an "algebraic view" of an operator. Algebraic
means, that the resulting operator acts on plain arrays,
instead of `FEFunctions`. This can be useful for solving
with external tools like `NLsolve.jl`.
See also [`FEOperator`](@ref).
"""
function get_algebraic_operator(feop::FEOperator)
  AlgebraicOpFromFEOp(feop)
end

struct AlgebraicOpFromFEOp <: NonlinearOperator
  feop::FEOperator
end

function allocate_residual(op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  allocate_residual(op.feop,u)
end

function residual!(b::AbstractVector,op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  residual!(b,op.feop,u)
end

function residual(op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  residual(op.feop,u)
end

function allocate_jacobian(op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  allocate_jacobian(op.feop,u)
end

function jacobian!(A::AbstractMatrix,op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  jacobian!(A,op.feop,u)
end

function jacobian(op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  jacobian(op.feop,u)
end

function residual_and_jacobian!(b::AbstractVector,A::AbstractMatrix,op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  residual_and_jacobian!(b,A,op.feop,u)
end

function residual_and_jacobian(op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  residual_and_jacobian(op.feop,u)
end

function zero_initial_guess(op::AlgebraicOpFromFEOp)
  trial = get_trial(op.feop)
  x = zero_free_values(trial)
end
