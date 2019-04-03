"""
Create 1-dim univariate polynomial basis of `UnivariateMonomialBasis` type
"""
function UnivariatePolynomialBasis(order::Int64)
  UnivariateMonomialBasis(order)
end

"""
Evaluate univariate monomial basis in a set of 1D points
"""
function (monomials::UnivariateMonomialBasis)(points::Vector{Float64})
  dbas = monomials.order+1
  c = Array{Float64,2}(undef, dbas, length(points))
  for (j,p) ∈ enumerate(points)
    for i= 1:dbas
      c[i,j] = p^(i-1)
    end
  end
  return c
end

"""
Function to be eliminated in the future. Compute the numder-th derivative of a monomial
at a set of 1D point
"""
function derivative(monomials::UnivariateMonomialBasis, numder::Int64,
                    points::Vector{Float64})
  dbas = monomials.order+1
  c = Array{Float64,2}(undef, dbas, length(points))
  for (j,p) ∈ enumerate(points)
    for i=1:size(c,1)
      c[i,j] = (i<=numder) ? 0.0 : prod([i-k-1 for k=0:numder-1])p^(i-numder-1)
    end
  end
  return c
end

"""
Provide a `TensorProductPolynomialBasis` for a vector `order` providing the order per
dimension
"""
function TensorProductPolynomialBasis(order::Vector{Int64})
  TensorProductPolynomialBasis([UnivariatePolynomialBasis(order[i]) for i=1:length(order)])
end

# @santiagobadia : If I put any comment to the function Documenter crashes
#"""
# Evaluate a `TensorProductPolynomialBasis` at an array of `Point{D}`
#"""
function (a::TensorProductPolynomialBasis)(points::Vector{Point{D}}) where {D}
  numdims = length(a.polynomials)
  @assert numdims == D
  c = Vector{Array{Float64,2}}(undef,D)
  for i ∈ 1:D
    pi = [p[i] for p ∈ points]
    c[i] = a.polynomials[i](pi)
  end
  # reinterpret(SVector{3,Float64}, x, (size(x,2),)) ?
  # @santiagobadia : Unfortunately, I have to do all this to get i-th coordinate
  orders = [(a.polynomials[i].order+1) for i ∈ 1:D]
  dims = tuple(orders...)
  A = Array{Float64,D}(undef,dims)
  B = Array{Float64,2}(undef,length(A),length(points))
  for j = 1:length(points)
    A.= 1
    tensorproduct!(A,c,j)
    B[:,j] = reshape(A, length(A))
  end
  return B
end

# function (a::TensorProductPolynomialBasis)(numders::Int64, x::Array{Float64,1})
#   numdims = length(a.polynomials)
#   @assert numdims == size(x,2) "Point dim and polynomial basis dim must be identical"
#   c = [[derivative(a.polynomials[i], j, x[:,i]) for j=1:numders] for i = 1:numdims]
#   orders = [(a.polynomials[i].order+1) for i = 1:numdims]
#   dims = tuple(orders...)
#   # A=ones(Float64,dims)
#   A = Array{Float64,length(c)}(undef,dims)
#   A.= 1.0
#   tensorproduct!(A,c)
#   A = reshape(A, length(A))
#   return A
# end

@generated function tensorproduct!(A::Array{Float64,N},c,ip=1) where {N}
  quote
    @nloops $N i A begin
    @nexprs $N j->( (@nref $N A i) *= c[j][i_j,ip])
  end
end
end

# @generated function tensorproductsquare!(A::Array{Array{Float64,N},M},c) where {N,M}
#   quote
#     @nloops $N i A begin
#     B = (@nref $N A i)
#     @nloops $M k B begin
#     @nexprs $N j ->((@nref $N B k) *= c[j][k_j,i_j])
#   end
# end
# end
# end

function gradient(a::TensorProductPolynomialBasis,
                  points::Vector{Point{D}}) where D
  @assert D == length(a.polynomials)
  x = reshape( reinterpret(Float64,points), (D,length(points)))
  c = [a.polynomials[i](x[i,:]) for i=1:D]
  dc = [derivative(a.polynomials[i], 1, x[i,:]) for i=1:D]
  orders = [(a.polynomials[i].order+1) for i=1:D]
  pldims = tuple(orders...)
  vals = Array{Float64,D}(undef,pldims)
  grad = Array{Float64,3}(undef,length(vals),size(x,2),D)
  for k=1:D
    d = [ (i==k) ? dc[i] : c[i] for i=1:D]
    for j = 1:size(x,2)
      vals .= 1.0
      tensorproduct!(vals,d,j)
      grad[:,j,k] = copy(reshape(vals,length(vals)))
    end
  end
  return grad
end
