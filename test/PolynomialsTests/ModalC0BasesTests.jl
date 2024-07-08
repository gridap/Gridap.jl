module ModalC0BasesTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
# using BenchmarkTools

import Gridap.Fields: Broadcasting

# Real-valued Q space with isotropic order

x1 = Point(0.0)
x2 = Point(0.5)
x3 = Point(1.0)

V = Float64
G = gradient_type(V,x1)
H = gradient_type(G,x1)

order = 3
a = fill(Point(-0.5),order+1)
b = fill(Point(2.5),order+1)
b1 = ModalC0Basis{1}(V,order,a,b)
∇b1 = Broadcasting(∇)(b1)
∇∇b1 = Broadcasting(∇)(∇b1)

@test evaluate(b1,[x1,x2,x3,]) ≈ [1.0 0.0 -0.0 0.0;
                                  0.5 0.5 -0.4330127018922193 0.18633899812498247;
                                  0.0 1.0 -0.0 -0.0]
@test evaluate(∇b1,[x1,x2,x3,]) ≈ G[(-1.0,) (1.0,) (-1.7320508075688772,) (1.4907119849998598,);
                                    (-1.0,) (1.0,) (-0.0,) (-0.37267799624996495,);
                                    (-1.0,) (1.0,) (1.7320508075688772,) (-0.0,)]
@test evaluate(∇∇b1,[x1,x2,x3,]) ≈ H[(0.0,) (0.0,) (3.4641016151377544,) (-5.962847939999439,);
                                     (0.0,) (0.0,) (3.4641016151377544,) (-1.4907119849998598,);
                                     (0.0,) (0.0,) (3.4641016151377544,) (2.9814239699997196,)]

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

end # module