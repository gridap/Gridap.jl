module PLambdaBasisTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using ForwardDiff
using StaticArrays
using BenchmarkTools

# TODO validation of precomputation in reference simplex

T = Float64

# 0D                                           0D #
D = 0
vertices = (Point{D,T}(),)
x = [vertices[1]]
x1 = x[1]

for k in 0:D
  for r in 1:4
    b = PmLambdaBasis(Val(D),T,r,k)
    evaluate(b,x)
    evaluate(Broadcasting(∇)(b),x)
    evaluate(Broadcasting(∇∇)(b),x)

    b = PmLambdaBasis(Val(D),T,r,k,vertices)
    evaluate(b,x)
    evaluate(Broadcasting(∇)(b),x)
    evaluate(Broadcasting(∇∇)(b),x)

    b = PLambdaBasis(Val(D),T,r,k)
    evaluate(b,x)
    evaluate(Broadcasting(∇)(b),x)
    evaluate(Broadcasting(∇∇)(b),x)

    b = PLambdaBasis(Val(D),T,r,k,vertices)
    evaluate(b,x)
    evaluate(Broadcasting(∇)(b),x)
    evaluate(Broadcasting(∇∇)(b),x)
  end
end

# 1D                                           1D #
D = 1
Pt = Point{D,T}
vertices = (Pt(.5),Pt(1.))
x = [xi for xi in vertices]
x1 = x[1]

r = 4
for k in 0:D
  b = PmLambdaBasis(Val(D),T,r,k)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  b = PmLambdaBasis(Val(D),T,r,k,vertices)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  b = PLambdaBasis(Val(D),T,r,k)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  b = PLambdaBasis(Val(D),T,r,k,vertices)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)
end


# 2D                                           2D #
D = 2
Pt = Point{D,T}
vertices = (Pt(0., 0.5),Pt(1.,0),Pt(.5,1.))
x = [xi for xi in vertices]
x1 = x[1]

r = 4
for k in 0:D
  b = PmLambdaBasis(Val(D),T,r,k)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  b = PmLambdaBasis(Val(D),T,r,k,vertices)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  b = PLambdaBasis(Val(D),T,r,k)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  b = PLambdaBasis(Val(D),T,r,k,vertices)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)
end

# 3D                                           3D #
D = 3
Pt = Point{D,T}
vertices = (Pt(0., 0., 0.5),Pt(1.,0,0),Pt(0,.5,0),Pt(0,.5,.5))
x = [xi for xi in vertices]
x1 = x[1]

r = 4
for k in 0:D
  b = PmLambdaBasis(Val(D),T,r,k)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  b = PmLambdaBasis(Val(D),T,r,k,vertices)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  b = PLambdaBasis(Val(D),T,r,k)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  b = PLambdaBasis(Val(D),T,r,k,vertices)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)
end


# 4D                                           4D #
D = 4
Pt = Point{D,T}
vertices = (Pt(0.,0.,0.,0.),Pt(0.,0.,0.,0.5),Pt(0.,1.,0,0),Pt(0.,0,.5,0),Pt(.5,1,1,1))
x = [xi for xi in vertices]
x1 = x[1]

r = 4
for k in (0,1,D-1,D)
  b = PmLambdaBasis(Val(D),T,r,k)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  b = PmLambdaBasis(Val(D),T,r,k,vertices)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  b = PLambdaBasis(Val(D),T,r,k)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  b = PLambdaBasis(Val(D),T,r,k,vertices)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)
end

end # module
