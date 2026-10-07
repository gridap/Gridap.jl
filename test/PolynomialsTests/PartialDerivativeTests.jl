module PartialDerivativeTests

using Test
using Gridap.Arrays
using Gridap.Fields
using Gridap.Polynomials
using Gridap.TensorValues

xs2 = [Point(0.3,-0.7), Point(0.1,0.4)]


###############################################
# Partial derivatives of a scalar valued basis
###############################################

# The gradient and hessian kernels fill the value type that `derivatives_type`
# prescribes, so the partials hold the same numbers as `FieldGradientArray`.

sb = MonomialBasis{2}(Float64, (2,2))

∂sb  = FieldPartialDerivativesArray{1}(sb)
∂²sb = FieldPartialDerivativesArray{2}(sb)

p1 = evaluate(∂sb, xs2)
p2 = evaluate(∂²sb, xs2)
g1 = evaluate(FieldGradientArray{1}(sb), xs2)
h2 = evaluate(FieldGradientArray{2}(sb), xs2)

@test eltype(p1) === derivatives_type(Float64, first(xs2), Val(1))
@test eltype(p2) === derivatives_type(Float64, first(xs2), Val(2))
@test size(p1) == size(g1) == size(p2)

@test all(p1[i,k][j] == g1[i,k][j] for i in axes(p1,1), k in axes(p1,2), j in 1:2)
# order 2 stores one component per non-decreasing derivation multi-index, and reads
# back the full symmetric hessian
@test all(p2[i,k][j1,j2] ≈ h2[i,k][j1,j2]
          for i in axes(p2,1), k in axes(p2,2), j1 in 1:2, j2 in 1:2)

# reaching the fused kernels is what keeps evaluation allocation free
c1 = return_cache(∂sb, xs2); evaluate!(c1, ∂sb, xs2)
c2 = return_cache(∂²sb, xs2); evaluate!(c2, ∂²sb, xs2)
@test iszero(@allocated evaluate!(c1, ∂sb, xs2))
@test iszero(@allocated evaluate!(c2, ∂²sb, xs2))

# `Broadcasting` transports the order into the array wrapper
@test Broadcasting(partial_derivatives)(sb) isa FieldPartialDerivativesArray{1}
@test Broadcasting(second_partial_derivatives)(sb) isa FieldPartialDerivativesArray{2}


###############################################
# Partial derivatives of a tensor valued basis
###############################################

# There `derivatives_type` is the product of the `SymNTensorValue` holding the
# derivatives of one basis polynomial with the value type of the basis. Its full
# tensor holds the same numbers as that of `FieldGradientArray`, which is dense, so
# comparing the two validates the dof order and the component layout at once.

# a component wise basis needs one row of orders per independent component of V
tensor_bases(V) = (
  MonomialBasis{2}(V,(2,2)),                                     # _cartprod_set_derivative!
  BernsteinBasis(Val(2),V,2),                                    # idem
  ModalC0Basis{2}(V,1),                                          # _set_derivative_mc0!
  CompWiseTensorPolyBasis{2}(Monomial,V,                         # _comp_wize_set_derivative!
                             fill(2,num_indep_components(V),2)))

for V in (ContraVectorValue{con2,Float64,2}, TensorValue{2,2,Float64},
          SymTensorValue{2,Float64}, SymTracelessTensorValue{2,Float64})

  for b in tensor_bases(V), Np in (1,2)
    ∂b = FieldPartialDerivativesArray{Np}(b)
    p  = evaluate(∂b, xs2)
    g  = evaluate(FieldGradientArray{Np}(b), xs2)

    @test eltype(p) === derivatives_type(V, first(xs2), Val(Np))
    @test size(p) == size(g)
    @test all(get_array(p[i,k]) == get_array(g[i,k]) for i in axes(p,1), k in axes(p,2))

    # reaching the fused kernels is what keeps evaluation allocation free
    cb = return_cache(∂b, xs2)
    evaluate!(cb, ∂b, xs2)
    @test iszero(@allocated evaluate!(cb, ∂b, xs2))

    # a single point is evaluated through the vector of points implementation
    cp = return_cache(∂b, first(xs2))
    @test evaluate!(cp, ∂b, first(xs2)) == p[1,:]
  end
end

# The Raviart-Thomas basis mixes a component wise block with dofs dense in V
rtb = RaviartThomasPolyBasis{2}(Monomial,Float64,1)
∂rtb = FieldPartialDerivativesArray{1}(rtb)
prt = evaluate(∂rtb, xs2)
grt = evaluate(FieldGradientArray{1}(rtb), xs2)
@test eltype(prt) === derivatives_type(VectorValue{2,Float64}, first(xs2), Val(1))
@test all(get_array(prt[i,k]) == get_array(grt[i,k]) for i in axes(prt,1), k in axes(prt,2))


###############################################
# Partial derivatives of a Barycentric PΛ basis
###############################################

# These bases build each dof as a sum of terms, every term the product of the
# derivatives of one scalar Bernstein polynomial with a constant form. The derivation
# indices are stored by the symmetric tensor of `derivatives_type`, so the product
# holds one component per non-decreasing derivation multi-index and per form component.
# A `k`=0 or `k`=`D` form is scalar valued, and its dofs then hold the derivation
# indices tensored with the order 0 `ScalarValue` carrying the value.

xs3 = [Point(0.3,-0.7,0.2), Point(0.1,0.4,0.3)]

bary_cases = (
  (BarycentricPmΛBasis(Val(2),Float64,2,0),              xs2), # scalar valued
  (BarycentricPΛBasis( Val(2),Float64,2,0),              xs2),
  (BarycentricPmΛBasis(Val(2),Float64,2,1),              xs2),
  (BarycentricPmΛBasis(Val(2),Float64,2,1;flavor=:BMM),  xs2), # divides the dofs by a multinomial
  (BarycentricPΛBasis( Val(2),Float64,2,1),              xs2),
  (BarycentricPmΛBasis(Val(2),Float64,2,2),              xs2), # scalar valued
  (BarycentricPmΛBasis(Val(3),Float64,2,1),              xs3),
  (BarycentricPmΛBasis(Val(3),Float64,2,2),              xs3),
  (BarycentricPΛBasis( Val(3),Float64,2,2),              xs3),
  (BarycentricPΛBasis( Val(3),Float64,2,3),              xs3)) # scalar valued

for (b,xs) in bary_cases, Np in (1,2)
  ∂b = FieldPartialDerivativesArray{Np}(b)

  if value_type(b) <: MultiValue
    @test_broken evaluate(∂b, xs) isa AbstractArray
    continue
  end

  p  = evaluate(∂b, xs)
  g  = evaluate(FieldGradientArray{Np}(b), xs)

  @test eltype(p) === derivatives_type(value_type(b), first(xs), Val(Np))
  @test size(p) == size(g)
  # the full tensor is that of the dense gradient, which a wrong component order or a
  # term built in the wrong type would break
  @test all(get_array(p[i,k]) == get_array(g[i,k]) for i in axes(p,1), k in axes(p,2))

  # reaching the fused kernels is what keeps evaluation allocation free
  cb = return_cache(∂b, xs)
  evaluate!(cb, ∂b, xs)
  @test iszero(@allocated evaluate!(cb, ∂b, xs))
end


end # module
