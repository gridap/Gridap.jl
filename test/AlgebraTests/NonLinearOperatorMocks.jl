module NonLinearOperatorMocks

export NonLinearOperatorMock

using Gridap
import Gridap: residual!, jacobian, jacobian!
import Gridap: create_in_domain

struct NonLinearOperatorMock <: NonLinearOperator end

function residual!(b::AbstractVector,::NonLinearOperatorMock,x::AbstractVector)
  b[1] = (x[1]-2)*(x[1]-1)
  b[2] = x[2]-3
end

function jacobian!(A::AbstractMatrix,::NonLinearOperatorMock,x::AbstractVector)
  A[1,1] = (x[1]-1) + (x[1]-2)
  A[1,2] = 0
  A[2,1] = 0
  A[2,2] = 1
end

function jacobian(op::NonLinearOperatorMock,x::AbstractVector)
  A = zeros(2,2)
  jacobian!(A,op,x)
  A
end

function create_in_domain(::NonLinearOperatorMock)
  zeros(2)
end

end # module
