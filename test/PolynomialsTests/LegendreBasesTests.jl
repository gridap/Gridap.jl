module LegendreBasisTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Fields: Broadcasting
using Gridap.Polynomials
using PolynomialBases
using ForwardDiff

@test isHierarchical(Legendre) == true

x1 = Point(0.0)
x2 = Point(0.5)
x3 = Point(1.0)
x = [x1, x2, x3]

V = Float64
G = gradient_type(V,x1)
H = gradient_type(G,x1)

function _legendre(n)
  t -> isone(n) ? 1. : sqrt(2*n-1)*jacobi(2t-1,n-1,0,0)
end
_∇(b) = t -> ForwardDiff.derivative(b, t)
_H(b) = t -> ForwardDiff.derivative(y -> ForwardDiff.derivative(b, y), t)

_bx_1D( order,x)   = [      _legendre(n)( xi[1])  for xi in x,  n in 1:order+1]
_Gbx_1D(order,x,G) = [ G(_∇(_legendre(n))(xi[1])) for xi in x,  n in 1:order+1]
_Hbx_1D(order,x,H) = [ H(_H(_legendre(n))(xi[1])) for xi in x,  n in 1:order+1]

order = 30
b = LegendreBasis(Val(1),V,order)

bx  = _bx_1D( order,x)
Gbx = _Gbx_1D(order,x,G)
Hbx = _Hbx_1D(order,x,H)

test_field_array(b,x,bx,≈, grad=Gbx, gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],≈, grad=Gbx[1,:], gradgrad=Hbx[1,:])

end # module
