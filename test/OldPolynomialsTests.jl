module PolynomialsTests
##
using Test
using Numa.FieldValues
using Numa.Polynomials
using Numa.Quadratures
using Numa.Polynomials: TensorProductPolynomialBasis
##

##

@testset "Mocks" begin

  include("PolynomialsTestsMocks.jl")

  basis = ShapeFunctionsScalarQua4()
  n = length(basis)
  @test n == 4


  quad = TensorProductQuadrature(orders=(2,2))
  points = coordinates(quad)
  values = Array{Float64,2}(undef, (n, length(points)) )
  evaluate!(basis,points,values)
  grad_basis = ∇*basis
  @test length(grad_basis) == 4

  grad_values = Array{VectorValue{2},2}(undef, (n, length(points)) )
  evaluate!(grad_basis,points,grad_values)

end

##

##
a = TensorProductPolynomialBasis([2,3])
p = Point{2}(2.0, 3.0)
b = a([p])
@test b ≈ [1.0 2.0 4.0 3.0 6.0 12.0 9.0 18.0 36.0 27.0 54.0 108.0]'
##

##
#Create a 1D monomial basis
monpol=UnivariatePolynomialBasis(2)
#Evaluate a monomial in a point
monpol([2.0])
#Evaluate derivatives of a monomial
a = [0.0]
using Numa.Polynomials: derivative
monpolder=derivative(monpol,2,[0.0])
@test monpolder[3] == 2.0
##

##
#Tensor product polynomial basis for multidimensional problems
multidpolb=TensorProductPolynomialBasis([2,3,2])
point = Point{3}(2.0,3.0,4.0)
mdpbval=multidpolb([point])
mdpbval[end]
@test (mdpbval[end]==(2.0^2*3.0^3*4.0^2))
##


##
#Evaluation of derivatives for monomials
polb = UnivariatePolynomialBasis(2)
#polb=TensorProductPolynomialBasis([2], basistype="Monomial")
polval=polb([3.0])
derval1=derivative(polb, 1, [3.0; 4.0])
derval2=derivative(polb, 2, [3.0; 4.0])
@test (derval1[end,1]==(2.0*3.0))
@test (derval2[end,1]==(2.0))
@test (derval1[end,2]==(2.0*4.0))
@test (derval2[end,2]==(2.0))
##

##
D = 2
x = [ Point{2}(2, 3), Point{2}(4, 5), Point{2}(7, 8), Point{2}(9, 10),
      Point{2}(11, 12) ]
typeof(x[1]) <: Point{2}
a=TensorProductPolynomialBasis([2,3])
B = a(x)
@test B[12,5] ≈ 209088
##

##
#Now multidimensional polynomial and gradients
points = [ Point{2}(3, 2), Point{2}(4, 6), Point{2}(5, 7) ]
polb=TensorProductPolynomialBasis([2,2])
grads = gradient(polb,points)
@test (grads[end,1,1]==24.0)
@test (grads[end,1,2]==36.0)
@test (grads[end,2,1]==288.0)
@test (grads[end,2,2]==192.0)
##

##
# Efficient evaluation of tensor product functions
orders=[2,3]
gps=(1,2)
quad = TensorProductQuadrature(orders=gps)
a=TensorProductPolynomialBasis(orders)
numdims = length(a.polynomials)
A = a(quad.coords)
@test size(A,1)== 12
@test A[4,1] ≈ -1/√3
@test A[4,2] ≈ 1/√3
##

# Evaluate tensor product-based computation of monomial basis gradients
##
orders=[2,3]
gps=(1,2)
quad = TensorProductQuadrature(orders=gps)
a=TensorProductPolynomialBasis(orders)
grad = gradient(a,quad.coords)
@test size(grad)==(12,2,2)
@test grad[7,1,2]==-1.1547005383792517
##

end # module PolynomialsTests
