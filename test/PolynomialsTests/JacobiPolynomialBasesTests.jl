module JacobiPolynomialBasisTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Fields: Broadcasting
using Gridap.Polynomials

# Real-valued Q space with isotropic order

x1 = Point(0.0)
x2 = Point(0.5)
x3 = Point(1.0)

V = Float64
G = gradient_type(V,x1)
H = gradient_type(G,x1)

order = 3
b1 = JacobiPolynomialBasis{1}(V,order)
∇b1 = Broadcasting(∇)(b1)
∇∇b1 = Broadcasting(∇)(∇b1)

@test evaluate(b1,[x1,x2,x3,]) ≈ [ 1.0 -1.7320508075688772 2.23606797749979 -2.6457513110645907;
                                   1.0 0.0 -1.118033988749895 -0.0;
                                   1.0 1.7320508075688772 2.23606797749979 2.6457513110645907 ]
@test evaluate(∇b1,[x1,x2,x3,]) ≈ G[ (0.0,) (3.4641016151377544,) (-13.416407864998739,) (31.74901573277509,); 
                                     (0.0,) (3.4641016151377544,) (0.0,) (-7.937253933193772,); 
                                     (0.0,) (3.4641016151377544,) (13.416407864998739,) (31.74901573277509,) ]
@test evaluate(∇∇b1,[x1,x2,x3,]) ≈ H[ (0.0,) (0.0,) (13.416407864998739,) (-79.37253933193772,); 
                                      (0.0,) (0.0,) (13.416407864998739,) (0.0,);
                                      (0.0,) (0.0,) (13.416407864998739,) (79.37253933193772,) ]

x1 = Point(0.0,0.0)
x2 = Point(0.5,0.5)
x3 = Point(1.0,1.0)
b2 = JacobiPolynomialBasis{2}(V,order)
∇b2 = Broadcasting(∇)(b2)
∇∇b2 = Broadcasting(∇)(∇b2)

G = gradient_type(V,x1)
H = gradient_type(G,x1)

@test evaluate(b2,[x1,x2,x3,]) ≈ [ 1.0 -1.7320508075688772 2.23606797749979 -2.6457513110645907 #=
                                =# -1.7320508075688772 2.9999999999999996 -3.872983346207417 #=
                                =# 4.58257569495584 2.23606797749979 -3.872983346207417 #=
                                =# 5.000000000000001 -5.916079783099617 -2.6457513110645907 #=
                                =# 4.58257569495584 -5.916079783099617 7.000000000000001; 
                                   1.0 0.0 -1.118033988749895 -0.0 0.0 0.0 -0.0 -0.0 #=
                                =# -1.118033988749895 -0.0 1.2500000000000002 0.0 -0.0 -0.0 0.0 0.0; 
                                   1.0 1.7320508075688772 2.23606797749979 2.6457513110645907 #= 
                                =# 1.7320508075688772 2.9999999999999996 3.872983346207417 #=
                                =# 4.58257569495584 2.23606797749979 3.872983346207417 #=
                                =# 5.000000000000001 5.916079783099617 2.6457513110645907 #= 
                                =# 4.58257569495584 5.916079783099617 7.000000000000001 ]
@test evaluate(∇b2,[x1,x2,x3,])[:,10] ≈ G[ (7.745966692414834, 23.2379000772445); 
                                           (-3.872983346207417, 0.0); 
                                           (7.745966692414834, 23.2379000772445) ]
@test evaluate(∇∇b2,[x1,x2,x3,])[:,10] ≈ H[ (0.0, -46.475800154489, -46.475800154489, -23.2379000772445); 
                                            (-0.0, 0.0, 0.0, 0.0);
                                            (0.0, 46.475800154489, 46.475800154489, 23.2379000772445) ]

end # module
