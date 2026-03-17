module ModalC0BasesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using StaticArrays
using ForwardDiff
using PolynomialBases
# using BenchmarkTools

import Gridap.Fields: Broadcasting

# Real-valued Q space with isotropic order

x1 = Point(0.0)
x2 = Point(0.5)
x3 = Point(1.0)
x = [x1, x2, x3]

V = Float64
G = gradient_type(V,x1)
H = gradient_type(G,x1)

function _modalC0(a,b,n)
    function(t)
      isone(n) && return  1. -t
      n==2     && return  t

      ξ = ( 2*t - ( a + b ) ) / ( b - a )
      return -sqrt(2*n-3)*t*(1-t)*jacobi(ξ,n-3,1,1)/(n-2)
    end
end
_∇(b) = t -> ForwardDiff.derivative(b, t)
_H(b) = t -> ForwardDiff.derivative(y -> ForwardDiff.derivative(b, y), t)


order = 30
a = fill(Point(-0.5),order+1)
b = fill(Point(2.5),order+1)
bm = ModalC0Basis{1}(V,order,a,b)

_bx_1D( order,x)   = [      _modalC0(a[1][1],b[1][1],n)( xi[1])  for xi in x,  n in 1:order+1]
_Gbx_1D(order,x,G) = [ G(_∇(_modalC0(a[1][1],b[1][1],n))(xi[1])) for xi in x,  n in 1:order+1]
_Hbx_1D(order,x,H) = [ H(_H(_modalC0(a[1][1],b[1][1],n))(xi[1])) for xi in x,  n in 1:order+1]

bx  = _bx_1D( order,x)
Gbx = _Gbx_1D(order,x,G)
Hbx = _Hbx_1D(order,x,H)

test_field_array(bm,x,bx,≈, grad=Gbx, gradgrad=Hbx)
test_field_array(bm,x[1],bx[1,:],≈, grad=Gbx[1,:], gradgrad=Hbx[1,:])

@test IndexStyle(bm) == IndexLinear()

# Validate generic 1D implem using CartProdPolyBasis

order = 3
a = fill(Point(0.),order+1)
b = fill(Point(1.),order+1)
b1 = ModalC0Basis{1}(V,order,a,b)
@test testvalue(typeof(b1)) isa typeof(b1)
b1u= CartProdPolyBasis(ModalC0,Val(1),V,order)

∇b1  = Broadcasting(∇)(b1)
∇b1u = Broadcasting(∇)(b1u)
∇∇b1 = Broadcasting(∇)(∇b1)
∇∇b1u= Broadcasting(∇)(∇b1u)

@test evaluate(b1,  [x1,x2,x3,]) ≈ evaluate(b1u,  [x1,x2,x3,])
@test evaluate(∇b1, [x1,x2,x3,]) ≈ evaluate(∇b1u, [x1,x2,x3,])
@test evaluate(∇∇b1,[x1,x2,x3,]) ≈ evaluate(∇∇b1u,[x1,x2,x3,])


x1 = Point(0.0,0.0)
x2 = Point(0.5,0.5)
x3 = Point(1.0,1.0)
a = [ Point(0.0,0.0),Point(0.0,0.0),Point(0.0,0.0),Point(0.0,0.0),
      Point(-0.5,2.5),Point(-0.5,2.5),Point(0.0,1.5),Point(0.0,1.5),
      Point(-1.0,-1.0),Point(-1.0,-1.0),Point(-1.0,-1.0),Point(-1.0,-1.0),
      Point(0.0,0.0),Point(0.0,0.0),Point(0.0,0.0),Point(0.0,0.0) ]
b = [ Point(1.0,1.0),Point(1.0,1.0),Point(1.0,1.0),Point(1.0,1.0),
      Point(2.5,2.5),Point(2.5,2.5),Point(1.5,1.5),Point(1.5,1.5),
      Point(-1.0,1.0),Point(-1.0,1.0),Point(-1.0,1.25),Point(-1.0,1.25),
      Point(1.0,1.0),Point(1.0,1.0),Point(1.0,1.0),Point(1.0,1.0) ]
b2 = ModalC0Basis{2}(V,order,a,b)
∇b2 = Broadcasting(∇)(b2)
∇∇b2 = Broadcasting(∇)(∇b2)

G = gradient_type(V,x1)
H = gradient_type(G,x1)

@test evaluate(b2,[x1,x2,x3,]) ≈ [ 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                                   0.25 0.25 0.25 0.25 -0.21650635094610965 0.09316949906249124 -0.21650635094610965 0.09316949906249124 -0.21650635094610965 -0.13975424859373686 -0.21650635094610965 -0.09316949906249122 0.18749999999999997 0.0 0.0 0.0;
                                   0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ]
@test evaluate(∇b2,[x1,x2,x3,])[:,10] ≈ G[ (0.0, 0.0,); (0.2795084971874737, -0.2795084971874737,); (0.0, 0.0,) ]
@test evaluate(∇∇b2,[x1,x2,x3,])[:,10] ≈ H[ (0.0, 0.0, 0.0, -4.47213595499958);
                                            (0.0, 0.5590169943749475, 0.5590169943749475, 1.118033988749895);
                                            (0.0, -2.23606797749979, -2.23606797749979, 0.0) ]

# Validate generic 2D implem using CartProdPolyBasis

order = 3
len_b2 = (order+1)^2
a = fill(Point(0.,0.), len_b2)
b = fill(Point(1.,1.), len_b2)

b2 = ModalC0Basis{2}(V,order,a,b)
b2u= CartProdPolyBasis(ModalC0,Val(2),V,order)
∇b2  = Broadcasting(∇)(b2)
∇b2u = Broadcasting(∇)(b2u)

b2x   = collect(eachcol(evaluate(b2,    [x1,x2,x3,])))
b2xu  = collect(eachcol(evaluate(b2u,   [x1,x2,x3,])))
∇b2x  = collect(eachcol(evaluate(∇b2,   [x1,x2,x3,])))
∇b2xu = collect(eachcol(evaluate(∇b2u,  [x1,x2,x3,])))

# re order basis polynomials as each basis has different ordering ...
b2x_perm   = b2x[  sortperm(b2x)[  invperm(sortperm(b2xu))]]
∇b2x_perm  = ∇b2x[ sortperm(∇b2x)[ invperm(sortperm(∇b2xu))]]

@test b2xu  == b2x_perm
@test ∇b2xu == ∇b2x_perm


# Misc

# Derivatives not implemented for symetric tensor types

D = 2
T = Float64
V = SymTensorValue{D,T}
G = gradient_type(V,x1)
s = MVector(0.,0.)
r = zeros(G, (1,1))
@test_throws ErrorException Polynomials._set_derivative_mc0!(r,1,s,0,0,V)

V = SymTracelessTensorValue{D,T}
G = gradient_type(V,x1)
r = zeros(G, (1,1))
@test_throws ErrorException Polynomials._set_derivative_mc0!(r,1,s,0,0,V)


end # module
