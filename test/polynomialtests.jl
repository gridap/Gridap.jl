
basis = ShapeFunctionsScalarQua4()

n = length(basis)
vt = valuetype(basis)

@test n == 4
@test vt == Float64

quad = TensorProductQuadrature{2}(orders=[2,2])

points = coordinates(quad)

values = Array{vt,2}(undef, (n, length(points)) )

evaluate!(basis,points,values)

grad_basis = gradient(basis)

grad_vt = valuetype(grad_basis)

@test grad_vt == VectorValue{2}

@test length(grad_basis) == 4

grad_values = Array{grad_vt,2}(undef, (n, length(points)) )

evaluate!(grad_basis,points,grad_values)


##
# a = TensorProductPolynomialBasis([1,2]; basistype="Lagrangian", nodestype="Equispaced")
# b = a([1.0 1.0])
# @test b ≈ [0.0 0.0 0.0 0.0 0.0 1.0]'
##
##
a = TensorProductPolynomialBasis([2,3])
p = Point{2}(2.0, 3.0)
b = a([p])
@test b ≈ [1.0 2.0 4.0 3.0 6.0 12.0 9.0 18.0 36.0 27.0 54.0 108.0]'
##
#Create a 1D monomial basis
monpol=UnivariatePolynomialBasis(2)
#Evaluate a monomial in a point
monpol([2.0])
#Evaluate derivatives of a monomial
monpolder=derivative(monpol,2,[0.0])

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
# x  = [2  3; 4 5; 7 8; 9 10; 11 12]
# a=TensorProductPolynomialBasis([2,3], basistype="Monomial")
# B = a(x)
# @test B[12,5] ≈ 209088
##
#Now multidimensional polynomial and gradients
x = [3.0 2.0; 4.0 6.0; 5.0 7.0]
polb=TensorProductPolynomialBasis([2,2], basistype="Monomial")
grads = gradient(polb,x)
@test (grads[end,1,1]==24.0)
@test (grads[end,1,2]==36.0)
@test (grads[end,2,1]==288.0)
@test (grads[end,2,2]==192.0)
##

# Efficient evaluation of tensor product functions
##
orders=[2,3]
gps=[1,2]
quad = TensorProductQuadrature{2}(orders=gps)
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
gps=[1,2]
quad = TensorProductQuadratureOld(gps)
a=TensorProductPolynomialBasis(orders, basistype="Monomial")
grad = gradient(a,quad.points)
@test size(grad)==(12,2,2)
@test grad[7,1,2]==-1.1547005383792517
##
