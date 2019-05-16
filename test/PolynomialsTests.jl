module PolynomialsTests

##
using Test
using Gridap
using Gridap.Polynomials
using Gridap.FieldValues

using Gridap.Polynomials: Basis
using Gridap.Polynomials: UnivariatePolynomialBasis
using Gridap.Polynomials: UnivariateMonomialBasis
using Gridap.Polynomials: evaluate
##

@testset "UnivariatePolynomialBasis" begin
  ##
  @test UnivariatePolynomialBasis <: Basis{1,ScalarValue}
  a = UnivariateMonomialBasis(2)
  @test typeof(a) <: Basis{1,ScalarValue}
  @test length(a) == 3
  #evaluate UnivariateMonomialBasis
  points = [Point{1}(1.0), Point{1}(2.0), Point{1}(3.0)]
  values = evaluate(a,points)
  @test values[3,3] == 9
  ##

  ##
  using Gridap.Polynomials: gradient
  grad = gradient(a)
  length(grad)
  using Gridap.Polynomials: evaluate
  gval = evaluate(grad, points)
  @test gval[3,3][1] == 6.0
  ##
end

@testset "TensorProductMonomialBasis" begin
  ##
  D=2
  orders=[2,3]
  tpmb = TensorProductMonomialBasis{D,ScalarValue}(orders)
  @test typeof(tpmb) <: Basis{2,ScalarValue}
  points = [ Point{2}(1.0, 1.0), Point{2}(1.0, 2.0), Point{2}(2.0, 2.0)]
  vals = evaluate(tpmb, points)
  @test vals[12,3] == 32.0
  ##

  ##
  orders=[1,1]
  tpmb = TensorProductMonomialBasis{D,ScalarValue}(orders)
  points = [ Point{2}(1.0, 1.0), Point{2}(1.0, 2.0)]
  T = ScalarValue
  this = tpmb
  points
  grbas = gradient(tpmb)
  vals = evaluate(grbas, points)
  vals
  @test vals[4,2][1] == 2.0
  ##

  ##
  orders=[1,1]
  tpmb = TensorProductMonomialBasis{D,VectorValue{D}}(orders)
  points = [ Point{2}(1.0, 1.0), Point{2}(1.0, 2.0)]
  T = VectorValue{D}
  this = tpmb
  points
  grbas = gradient(tpmb)
  v = evaluate(grbas, points)
  @test v[8,2][1,2] == 2.0
  ##

  ##
  orders=[2,3]
  tpmb = TensorProductMonomialBasis{D,VectorValue{D}}(orders)
  @test typeof(tpmb) <: Basis{2,VectorValue{D}}
  points = [ Point{2}(1.0, 1.0), Point{2}(1.0, 2.0), Point{2}(2.0, 2.0)]
  vals = evaluate(tpmb, points)
  @test vals[12,3][1] == 32.0
  @test vals[12,3][2] == 0.0
  @test vals[24,3][1] == 0.0
  @test vals[24,3][2] == 32.0
  ##

  ##
  orders=[2,3]
  tpmb = TensorProductMonomialBasis{D,TensorValue{D,D*D}}(orders)
  @test typeof(tpmb) <: Basis{2,TensorValue{D,D*D}}
  points = [ Point{2}(1.0, 1.0), Point{2}(1.0, 2.0), Point{2}(2.0, 2.0)]
  vals = evaluate(tpmb, points)
  @test vals[12,3][1,1] == 32.0
  @test vals[12,3][2,1] == 0.0
  @test vals[12,3][1,2] == 0.0
  @test vals[12,3][2,2] == 0.0
  @test vals[24,3][1,1] == 0.0
  @test vals[24,3][2,1] == 32.0
  @test vals[24,3][1,2] == 0.0
  @test vals[24,3][2,2] == 0.0
  @test vals[36,3][1,1] == 0.0
  @test vals[36,3][2,1] == 0.0
  @test vals[36,3][1,2] == 32.0
  @test vals[36,3][2,2] == 0.0
  @test vals[48,3][1,1] == 0.0
  @test vals[48,3][2,1] == 0.0
  @test vals[48,3][1,2] == 0.0
  @test vals[48,3][2,2] == 32.0
  ##
end

end #module PolynomialsTests
