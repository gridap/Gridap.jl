
"""
    abstract type LinearSolver <: GridapType end

- [`symbolic_setup(::LinearSolver,mat::AbstractMatrix)`](@ref)
- [`test_linear_solver`](@ref)

"""
abstract type LinearSolver <: GridapType end

"""
    symbolic_setup(::LinearSolver,mat::AbstractMatrix) -> SymbolicSetup
"""
function symbolic_setup(::LinearSolver,mat::AbstractMatrix)
  @abstractmethod
end

"""
    solve(ls::LinearSolver,A::AbstractMatrix,b::AbstractVector)
"""
function solve(ls::LinearSolver,A::AbstractMatrix,b::AbstractVector)
  ss = symbolic_setup(ls,A)
  ns = numerical_setup(ss,A)
  x = similar(b)
  solve!(x,ns,b)
  x
end

"""
    abstract type SymbolicSetup <: GridapType end

- [`numerical_setup(::SymbolicSetup,mat::AbstractMatrix)`](@ref)
"""
abstract type SymbolicSetup <: GridapType end

"""
    numerical_setup(::SymbolicSetup,mat::AbstractMatrix) -> NumericalSetup
"""
function numerical_setup(::SymbolicSetup,mat::AbstractMatrix)
  @abstractmethod
end

"""
    abstract type NumericalSetup <: GridapType end

- [`numerical_setup!(::NumericalSetup,mat::AbstractMatrix)`](@ref)
- [`solve!(x::AbstractVector,::NumericalSetup,b::AbstractVector)`](@ref)
"""
abstract type NumericalSetup <: GridapType end

"""
    numerical_setup!(::NumericalSetup,mat::AbstractMatrix)
"""
function numerical_setup!(::NumericalSetup,mat::AbstractMatrix)
  @abstractmethod
end

"""
    solve!(x::AbstractVector,::NumericalSetup,b::AbstractVector)
"""
function solve!(x::AbstractVector,::NumericalSetup,b::AbstractVector)
  @abstractmethod
end

"""
    test_linear_solver(
      ls::LinearSolver,
      A::AbstractMatrix,
      b::AbstractVector,
      x::AbstractVector)
"""
function test_linear_solver(
  ls::LinearSolver,
  A::AbstractMatrix,
  b::AbstractVector,
  x::AbstractVector)

  y = solve(ls,A,b)
  @test x ≈ y

  ss = symbolic_setup(ls,A)
  ns = numerical_setup(ss,A)
  numerical_setup!(ns,A)
  solve!(y,ns,b)
  @test x ≈ y

end

# Concrete implementations

"""
    struct LUSolver <: LinearSolver end

Wrapper of the LU solver available in julia
"""
struct LUSolver <: LinearSolver end

struct LUSymbolicSetup <: SymbolicSetup end

mutable struct LUNumericalSetup{F} <: NumericalSetup
  factors::F
end

symbolic_setup(::LUSolver,mat::AbstractMatrix) = LUSymbolicSetup()

numerical_setup(::LUSymbolicSetup,mat::AbstractMatrix) = LUNumericalSetup(lu(mat))

function numerical_setup!(ns::LUNumericalSetup, mat::AbstractMatrix)
  fac = lu(mat)
  ns.factors = fac
end

function solve!(
  x::AbstractVector,ns::LUNumericalSetup,b::AbstractVector)
  y = ns.factors\b # the allocation of y can be avoided
  x .= y
end

"""
    struct BackslashSolver <: LinearSolver end

Wrapper of the backslash solver available in julia
This is typically faster than LU for a single solve
"""
struct BackslashSolver <: LinearSolver end

struct BackslashSymbolicSetup <: SymbolicSetup end

mutable struct BackslashNumericalSetup{T<:AbstractMatrix} <: NumericalSetup
  A::T
end

symbolic_setup(::BackslashSolver,mat::AbstractMatrix) = BackslashSymbolicSetup()

numerical_setup(::BackslashSymbolicSetup,mat::AbstractMatrix) = BackslashNumericalSetup(mat)

function numerical_setup!(ns::BackslashNumericalSetup, mat::AbstractMatrix)
  ns.A = mat
end

function solve!(
  x::AbstractVector,ns::BackslashNumericalSetup,b::AbstractVector)
   copyto!(x, ns.A\b)
end
