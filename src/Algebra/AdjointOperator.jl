
"""
struct AffineFEOperatorWithAdjoint{A<:AbstractMatrix,B<:AbstractVector} <: NonlinearOperator
  matrix::A
  matrix_adj::A
  vector::B
  vector_adj::B
end
"""
struct AffineOperatorWithAdjoint{A<:AbstractMatrix,B<:AbstractVector} <: NonlinearOperator
  matrix::A
  matrix_adj::A
  vector::B
  vector_adj::B
end

"""
get_matrix(operator)

Return the matrix corresponding to the assembled left hand side of the operator.
This matrix incorporates all boundary conditions and constraints.
"""
get_matrix(op::AffineOperatorWithAdjoint) = op.matrix
"""
get_adjoint_matrix(operator)

Return the matrix corresponding to the assembled left hand side of the adjoint operator.
This matrix incorporates all boundary conditions and constraints.
"""
get_adjoint_matrix(op::AffineOperatorWithAdjoint) = op.matrix_adj

"""
get_vector(operator)

Return the vector corresponding to the assembled right hand side of the operator.
This vector includes all boundary conditions and constraints.
"""
get_vector(op::AffineOperatorWithAdjoint) = op.vector

"""
get_adjoint_vector(operator)

Return the vector corresponding to the assembled right hand side of the adjoint operator.
This vector includes all boundary conditions and constraints.
"""
get_adjoint_vector(op::AffineOperatorWithAdjoint) = op.vector_adj

function residual!(
  b::Tuple{AbstractVector,AbstractVector},
  op::AffineOperatorWithAdjoint,
  x::Tuple{AbstractVector,AbstractVector})
  mul!(b[1],op.matrix,x[1])
  mul!(b[2],op.matrix_adj,x[2])
  add_entries!(b[1],op.vector,-)
  add_entries!(b[2],op.vector_adj,-)
  b
end

function jacobian!(
  A::Tuple{AbstractMatrix,AbstractMatrix},
  op::AffineOperatorWithAdjoint,
  x::Tuple{AbstractVector,AbstractVector})
  copy_entries!(A[1],op.matrix)
  copy_entries!(A[2],op.matrix_adj)
  A
end

function jacobian(
  op::AffineOperatorWithAdjoint,
  x::Tuple{AbstractVector,AbstractVector})
  (op.matrix, op.matrix_adj)
end

function zero_initial_guess(op::AffineOperatorWithAdjoint)
  x = allocate_in_domain(typeof(op.vector),op.matrix)
  x_adj = allocate_in_domain(typeof(op.vector_adj),op.matrix_adj)
  fill_entries!(x,zero(eltype(x)))
  fill_entries!(x_adj,zero(eltype(x_adj)))
  (x, x_adj)
end

function allocate_residual(
  op::AffineOperatorWithAdjoint,
  x::Tuple{AbstractVector,AbstractVector})
  r = allocate_in_range(typeof(op.vector),op.matrix)
  r_adj = allocate_in_range(typeof(op.vector_adj),op.matrix_adj)
  (r,r_adj)
end

function allocate_jacobian(
  op::AffineOperatorWithAdjoint,
  x::Tuple{AbstractVector,AbstractVector})
  (op.matrix, op.matrix_adj)
end

"""
solve(ls::LinearSolver,A::Tuple{AbstractMatrix,AbstractMatrix},b::Tuple{AbstractVector,AbstractVector})
"""
function solve(
  ls::LinearSolver,
  A::Tuple{AbstractMatrix,AbstractMatrix},
  b::Tuple{AbstractVector,AbstractVector})
  ss = symbolic_setup(ls,A)
  ns = numerical_setup(ss,A)
  x = similar(b)
  solve!(x,ns,b)
  x
end

function solve!(x::Tuple{AbstractVector,AbstractVector},
  ls::LinearSolver,
  op::NonlinearOperator,
  cache::Nothing)
  fill_entries!(x,zero(eltype(x)))
  b = residual(op, x)
  A = jacobian(op, x)
  ss = symbolic_setup(ls, A)
  ns = numerical_setup(ss,A)
  scale_entries!(b,-1)
  solve!(x,ns,b)
  LinearSolverCache(A,b,ns)
end

function solve!(x::Tuple{AbstractVector,AbstractVector},
  ls::LinearSolver,
  op::NonlinearOperator,
  cache)
  fill_entries!(x,zero(eltype(x)))
  b = cache.b
  A = cache.A
  ns = cache.ns
  residual!(b, op, x)
  numerical_setup!(ns,A)
  scale_entries!(b,-1)
  solve!(x,ns,b)
  cache
end

function solve!(
  x::Tuple{AbstractVector,AbstractVector},
  ls::LinearSolver,
  op::AffineOperatorWithAdjoint,
  cache::Nothing)
  A = (op.matrix, op.matrix_adj)
  b = (op.vector, op.vector_adj)
  ss = symbolic_setup(ls,A)
  ns = numerical_setup(ss,A)
  solve!(x,ns,b)
  ns
end

function solve!(
  x::Tuple{AbstractVector,AbstractVector},
  ls::LinearSolver,
  op::AffineOperatorWithAdjoint,
  cache)
  newmatrix = true
  solve!(x,ls,op,cache,newmatrix)
  cache
end

"""
solve!(
x::Tuple{AbstractVector,AbstractVector},
ls::LinearSolver,
op::AffineOperatorWithAdjoint,
cache,
newmatrix::Bool)
"""
function solve!(
  x::Tuple{AbstractVector,AbstractVector},
  ls::LinearSolver,
  op::AffineOperatorWithAdjoint,
  cache,
  newmatrix::Bool)
  A = (op.matrix, op.matrix_adj)
  b = (op.vector, op.vector_adj)
  ns = cache
  if newmatrix
    numerical_setup!(ns,A)
  end
  solve!(x,ns,b)
  cache
end

function solve!(
  x::Tuple{AbstractVector,AbstractVector},
  ls::LinearSolver,
  op::AffineOperatorWithAdjoint,
  cache::Nothing,
  newmatrix::Bool)
  solve!(x,ls,op,cache)
end

# Concrete implementations

struct LUSymbolicSetupWithAdjoint <: SymbolicSetup end

mutable struct LUNumericalSetupWithAdjoint{F} <: NumericalSetup
  factors::F
  factors_adj
end

symbolic_setup(::LUSolver,mat::Tuple{AbstractMatrix,AbstractMatrix}) = LUSymbolicSetupWithAdjoint()

function numerical_setup(
  ::LUSymbolicSetupWithAdjoint,
  mat::Tuple{AbstractMatrix,AbstractMatrix})
  fac = lu(mat[1])
  fac_adj = transpose(fac)
  LUNumericalSetupWithAdjoint(fac,fac_adj)
end

function numerical_setup!(
  ns::LUNumericalSetupWithAdjoint,
  mat::Tuple{AbstractMatrix,AbstractMatrix})
  fac = lu(mat[1])
  fac_adj = transpose(fac)
  ns.factors = fac
  ns.factors_adj = fac_adj
end

function solve!(
  x::Tuple{AbstractVector,AbstractVector},
  ns::LUNumericalSetupWithAdjoint,
  b::Tuple{AbstractVector,AbstractVector})
  y = ns.factors\b[1] # the allocation of y can be avoided
  y_adj = ns.factors_adj\b[2]
  copy_entries!(x[1],y)
  copy_entries!(x[2],y_adj)
  x
end

struct BackslashSymbolicSetupWithAdjoint <: SymbolicSetup end

mutable struct BackslashNumericalSetupWithAdjoint{T<:AbstractMatrix} <: NumericalSetup
  A::T
  A_adj::T
end

symbolic_setup(::BackslashSolver,mat::Tuple{AbstractMatrix,AbstractMatrix}) = BackslashSymbolicSetupWithAdjoint()

function numerical_setup(
  ::BackslashSymbolicSetup,
  mat::Tuple{AbstractMatrix,AbstractMatrix})
  BackslashNumericalSetup(mat[1],mat[2])
end

function numerical_setup!(
  ns::BackslashNumericalSetup,
  mat::Tuple{AbstractMatrix,AbstractMatrix})
  ns.A = mat[1]
  ns.A_adj = mat[2]
end

function solve!(
  x::Tuple{AbstractVector,AbstractVector},
  ns::BackslashNumericalSetup,
  b::Tuple{AbstractVector,AbstractVector})
  copy_entries!(x, (ns.A\b[1],ns.A_adj\b[2])) # Not sure how to optimize the "\" solver
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
  A::Tuple{AbstractMatrix,AbstractMatrix},
  b::Tuple{AbstractVector,AbstractVector},
  x::Tuple{AbstractVector,AbstractVector})

  y = solve(ls,A,b)
  @test x[1] ≈ y[1]
  @test x[2] ≈ y[2]

  ss = symbolic_setup(ls,A)
  ns = numerical_setup(ss,A)
  numerical_setup!(ns,A)
  solve!(y,ns,b)
  @test x[1] ≈ y[1]
  @test x[2] ≈ y[2]

  op = AffineOperatorWithAdjoint(A[1],A[2],b[1],b[2])
  x0 = copy(x)
  test_nonlinear_solver(ls,op,x0,x)

end
