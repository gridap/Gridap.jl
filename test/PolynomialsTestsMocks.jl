
using Numa.Helpers
using Numa.FieldValues
using Numa.Polynomials

import Numa.Polynomials: evaluate!, gradient

struct ShapeFunctionsScalarQua4 <: MultivariatePolynomialBasis{2,Float64} end

struct GradShapeFunctionsScalarQua4 <: MultivariatePolynomialBasis{2,VectorValue{2}} end

Base.length(::ShapeFunctionsScalarQua4) = 4

function evaluate!(
  ::ShapeFunctionsScalarQua4,points::Array{Point{2},1},v::Array{Float64,2})
  for (i,point) in enumerate(points)
    xi = point[1]
    eta = point[2]
    v[1,i] = (1-xi)*(1-eta)/4.0
    v[2,i] = (1+xi)*(1-eta)/4.0
    v[3,i] = (1-xi)*(1+eta)/4.0
    v[4,i] = (1+xi)*(1+eta)/4.0
  end
end

gradient(::ShapeFunctionsScalarQua4) = GradShapeFunctionsScalarQua4()

Base.length(::GradShapeFunctionsScalarQua4) = 4

function evaluate!(
  ::GradShapeFunctionsScalarQua4,points::Array{Point{2},1},v::Array{VectorValue{2},2})
  for (i,point) in enumerate(points)
    xi = point[1]
    eta = point[2]
    v[1,i] = VectorValue{2}( (eta-1)/4.0, (xi-1)/4.0 )
    v[2,i] = VectorValue{2}( (1-eta)/4.0,-(1+xi)/4.0 )
    v[3,i] = VectorValue{2}(-(1+eta)/4.0, (1-xi)/4.0 )
    v[4,i] = VectorValue{2}( (1+eta)/4.0, (1+xi)/4.0 )
  end
end

gradient(::GradShapeFunctionsScalarQua4) = @notimplemented

