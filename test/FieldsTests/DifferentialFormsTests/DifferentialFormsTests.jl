module DifferentialFormsTests

using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: FieldPartialDerivatives, FieldPartialDerivativesArray
using Gridap.Fields: OperationField
using FillArrays: Fill
using Test
import LinearAlgebra

# ω1 = 3x¹ dx¹ + 9(x¹)² dx², ω2 swaps the two components
ω1fun(x) = CoVectorValue(3.0*x[1], 9.0*x[1]^2)
ω2fun(x) = CoVectorValue(9.0*x[1]^2, 3.0*x[1])
ω1 = GenericField(ω1fun)
ω2 = GenericField(ω2fun)

# ── Single-point ──────────────────────────────────────────────────────────────

x = Point(1.0,1.0)

cache = return_cache(ω1,x)
ω1x = evaluate!(cache,ω1,x)
@test ω1x.data == (3.0,9.0)

cache = return_cache(ω2,x)
ω2x = evaluate!(cache,ω2,x)
@test ω2x.data == (9.0,3.0)

k = ω1 ∧ ω2
kx = evaluate(k,x)
@test kx == ω1x ∧ ω2x
@test kx.data == (-72.0,)

k = ω1 ∧ ω1
kx = evaluate(k,x)
@test kx == ω1x ∧ ω1x
@test kx.data == (0.0,)

# ── Multi-point evaluation ────────────────────────────────────────────────────

xs = [Point(1.0,1.0), Point(2.0,1.0), Point(1.0,2.0)]

cache_v = return_cache(ω1, xs)
vals    = evaluate!(cache_v, ω1, xs)
@test vals isa AbstractVector
@test length(vals) == 3
@test vals[1].data == (3.0, 9.0)
@test vals[2].data == (6.0, 36.0)  # At (2,1): 3x¹ = 6, 9(x¹)² = 36
@test vals[3].data == (3.0, 9.0)   # At (1,2): the components ignore x²

# cache reuse with a different number of points
xs2   = [Point(0.5,0.5)]
vals2 = evaluate!(cache_v, ω1, xs2)
@test length(vals2) == 1
@test vals2[1].data == (1.5, 2.25)

# ── FieldPartialDerivatives ───────────────────────────────────────────────────
# The wrapper defers to the same numbers as the `Function` method, which nests
# `ForwardDiff.jacobian` directly.

∂ω1 = partial_derivatives(ω1)
@test ∂ω1 isa FieldPartialDerivatives{1}
@test ∂ω1 === partial_derivatives(ω1)

ref1 = partial_derivatives(ω1fun, x)
@test evaluate(∂ω1, x) === ref1
@test (@inferred evaluate(∂ω1, x)) === ref1
@test typeof(evaluate(∂ω1, x)) === derivatives_type(typeof(ω1x), x, Val(1))

∂²ω1 = second_partial_derivatives(ω1)
@test ∂²ω1 isa FieldPartialDerivatives{2}
@test evaluate(∂²ω1, x) === second_partial_derivatives(ω1fun, x)

# Multi point evaluation goes through the generic array method
@test evaluate(∂ω1, xs) == [evaluate(∂ω1, xi) for xi in xs]

# Order 0 is the value itself, promoted by the coordinate eltype
@test evaluate(PartialDerivatives{0}()(ω1), x) === ω1x*one(eltype(x))

# Derivation orders add up, the mixed partials of order N+M being held by one
# fully symmetric derivation tensor
@test partial_derivatives(∂ω1) === ∂²ω1
@test evaluate(partial_derivatives(∂ω1), x) === evaluate(∂²ω1, x)

# Order 0 leaves the partials alone, it does not promote their eltype like it
# does the value of an undifferentiated field
@test PartialDerivatives{0}()(∂ω1) === ∂ω1

# An array of fields is differentiated through `Broadcasting`
@test_throws ErrorException partial_derivatives([ω1, ω1])

# A product with a constant field keeps only the term that differentiates the
# other factor
opf = Operation(*)(ConstantField(2.0), ω1)
@test opf isa OperationField
@test evaluate(partial_derivatives(opf), x) ===
      partial_derivatives(y -> 2.0*ω1fun(y), x)

# Partials of a constant vanish without differentiating anything
cf = ConstantField(ω1x)
@test iszero(evaluate(partial_derivatives(cf), x))
@test evaluate(partial_derivatives(cf), x) ===
      zero(derivatives_type(typeof(ω1x), x, Val(1)))
@test evaluate(PartialDerivatives{0}()(cf), x) === ω1x

# Derivation is linear, so the partials of a sum are the sum of the partials
for op in (+,-)
  @test evaluate(partial_derivatives(Operation(op)(ω1,ω2)), x) ===
        partial_derivatives(y -> op(ω1fun(y),ω2fun(y)), x)
end
@test evaluate(partial_derivatives(ω1+ω2), x) ===
      partial_derivatives(y -> ω1fun(y)+ω2fun(y), x)

# A linear operation with constant coefficients commutes with derivation at every
# order, so it acts on the value factor of the partials alone
for Np in (1,2)
  @test evaluate(PartialDerivatives{Np}()(Operation(flat_hodge_star)(ω1)), x) ===
        PartialDerivatives{Np}()(y -> flat_hodge_star(ω1fun(y)), x)
end

# A product rule is a two term sum, each term acting on the value factor of one
# argument's partials
sfun(y) = 2.0*y[1]^2
ufun(y) = ContraVectorValue(y[1], 2.0*y[2])
wfun(y) = CoVectorValue(3.0, y[1])
for (op,f1fun,f2fun) in (
  (*,sfun,ω1fun), (⋅,ufun,wfun), (⊗,ufun,wfun), (∧,ω1fun,ω2fun)
)
  opf2 = Operation(op)(GenericField(f1fun), GenericField(f2fun))
  @test evaluate(partial_derivatives(opf2), x) ===
        partial_derivatives(y -> op(f1fun(y),f2fun(y)), x)
end

@test_broken evaluate(partial_derivatives(Operation(⊙)(GenericField(ufun), GenericField(wfun))), x) ===
             partial_derivatives(y -> ufun(y)⊙wfun(y), x)

# A contraction of the two factors is scalar valued, and its partials keep the
# product shape with a `ScalarValue` second factor
@test typeof(evaluate(partial_derivatives(Operation(⋅)(GenericField(ufun),GenericField(wfun))), x)) ===
      derivatives_type(typeof(ufun(x)⋅wfun(x)), x, Val(1))

# The order 2 rule adds the mixed term pairing the two first derivatives, which
# is the only non zero one when both factors are affine
for (op,f1fun,f2fun) in (
  (*,sfun,ω1fun), (⋅,ufun,wfun), (⊗,ufun,wfun), (∧,ω1fun,ω2fun)
)
  opf3 = Operation(op)(GenericField(f1fun), GenericField(f2fun))
  @test evaluate(second_partial_derivatives(opf3), x) ===
        second_partial_derivatives(y -> op(f1fun(y),f2fun(y)), x)
end

@test_broken evaluate(second_partial_derivatives(Operation(⊙)(GenericField(ufun), GenericField(wfun))), x) ===
             second_partial_derivatives(y -> ufun(y)⊙wfun(y), x)

# The two derivation indices share one fully symmetric tensor, which is not full
# rank at order 2, so the mixed term has to build that container itself
@test typeof(evaluate(second_partial_derivatives(Operation(⋅)(GenericField(ufun),GenericField(wfun))), x)) ===
      derivatives_type(typeof(ufun(x)⋅wfun(x)), x, Val(2))

# In 3D the two derivation indices have six independent pairs, so the order they
# are stored in is no longer trivial
x3 = Point(0.3,0.5,0.9)
s3fun(y) = 2.0*y[1]^2 + y[2]*y[3]
f3fun(y) = CoVectorValue(y[1]*y[2], 3.0*y[3]^2, y[1]+y[2])
g3fun(y) = CoVectorValue(y[2]^2, y[1], 4.0)
u3fun(y) = ContraVectorValue(y[1]^2, y[2]*y[3], y[1])
w3fun(y) = CoVectorValue(y[3], y[1]*y[2], 2.0)
for (op,f1fun,f2fun) in ((⋅,u3fun,w3fun), (⊗,u3fun,w3fun), (∧,f3fun,g3fun))
  opf4 = Operation(op)(GenericField(f1fun), GenericField(f2fun))
  @test evaluate(second_partial_derivatives(opf4), x3) ===
        second_partial_derivatives(y -> op(f1fun(y),f2fun(y)), x3)
end

# Both sides sum the scalar factor's terms in their own order, so one component
# of this one differs in its last bit
opf5 = Operation(*)(GenericField(s3fun), GenericField(f3fun))
@test evaluate(second_partial_derivatives(opf5), x3) ≈
      second_partial_derivatives(y -> s3fun(y)*f3fun(y), x3)
@test typeof(evaluate(second_partial_derivatives(opf5), x3)) ===
      derivatives_type(typeof(s3fun(x3)*f3fun(x3)), x3, Val(2))

# Above order 2 the rule is a multi index Leibniz sum, and a product rule needs
# exactly two factors
@test_throws ErrorException PartialDerivatives{3}()(Operation(*)(ConstantField(2.0), ω1))
@test_throws ErrorException partial_derivatives(Operation(*)(ConstantField(2.0), ω1, ω2))

# A scalar factor's partials have a `ScalarValue` value factor, which scales an
# order 0 form without dropping its embedding dimension
z0fun(y) = ExteriorFormValue{0,2}((2.0*y[1],))
s0fun(y) = 3.0*y[2]
opf0 = Operation(*)(GenericField(s0fun), GenericField(z0fun))
@test evaluate(partial_derivatives(opf0), x) ===
      partial_derivatives(y -> s0fun(y)*z0fun(y), x)
@test typeof(evaluate(partial_derivatives(opf0), x)) ===
      derivatives_type(typeof(s0fun(x)*z0fun(x)), x, Val(1))

# Zero and void fields carry the derivation through their wrapper
@test evaluate(partial_derivatives(ZeroField(ω1)), x) ===
      zero(derivatives_type(typeof(ω1x), x, Val(1)))
@test evaluate(partial_derivatives(Fields.VoidField(ω1,false)), x) ===
      evaluate(partial_derivatives(ω1), x)

# The partials of an affine field are its constant gradient, and it has no curvature
h = AffineField(TensorValue(1.0,2.0,3.0,4.0), Point(1.0,1.0))
@test evaluate(partial_derivatives(h), x) ===
      partial_derivatives(y -> y⋅h.gradient + h.origin, x)
@test iszero(evaluate(second_partial_derivatives(h), x))

# ── FieldPartialDerivativesArray ──────────────────────────────────────────────

fa  = Fill(ω1, 3)
∂fa = Broadcasting(partial_derivatives)(fa)
@test ∂fa isa FieldPartialDerivativesArray{1}
@test size(∂fa) == size(fa)
@test axes(∂fa) == axes(fa)
@test IndexStyle(typeof(∂fa)) === IndexStyle(typeof(fa))
@test eltype(∂fa) === typeof(∂ω1)
@test ∂fa[2] === partial_derivatives(fa[2])
@test evaluate(∂fa[2], x) === ref1
@test typeof(testvalue(typeof(∂fa))) === typeof(∂fa)

# The order is a type parameter of the operator, so `Broadcasting` carries it into
# the array rather than being pinned to the first derivatives
∂²fa = Broadcasting(second_partial_derivatives)(fa)
@test ∂²fa isa FieldPartialDerivativesArray{2}
@test ∂²fa[2] === second_partial_derivatives(fa[2])
@test evaluate(∂²fa[2], x) === second_partial_derivatives(ω1fun, x)
@test typeof(testvalue(typeof(∂²fa))) === typeof(∂²fa)
@test lazy_map(Broadcasting(second_partial_derivatives), fa)[1] === second_partial_derivatives(fa[1])

# Orders add up at the array level too
∂∂fa = Broadcasting(partial_derivatives)(∂fa)
@test ∂∂fa isa FieldPartialDerivativesArray{2}
@test ∂∂fa === ∂²fa

# ── exterior_derivative ───────────────────────────────────────────────────────
# dω1 = (∂₁(9(x¹)²) - ∂₂(3x¹)) dx¹∧dx² = 18x¹ dx¹∧dx²

dω1 = exterior_derivative(ω1)
@test evaluate(dω1, x) === partials_to_exteriorderivative(ref1)
@test evaluate(dω1, x) isa ExteriorFormValue{2,cov2}
@test evaluate(dω1, x).data[1] ≈ 18.0   atol=1e-12

dvals = evaluate(dω1, xs)
@test length(dvals) == 3
@test dvals[1].data[1] ≈ 18.0   atol=1e-12   # x¹ = 1
@test dvals[2].data[1] ≈ 36.0   atol=1e-12   # x¹ = 2
@test dvals[3].data[1] ≈ 18.0   atol=1e-12   # x¹ = 1

@test evaluate(∇∧ω1, x) === evaluate(dω1, x)  # ∇∧ spells the same operator
@test evaluate(𝑑(ω1), x) === evaluate(dω1, x)

# d of a 0-form proxied by a Number
sfun(x) = sin(x[1])*exp(x[2])
sfld = GenericField(sfun)
@test evaluate(exterior_derivative(sfld), x) ===
      partials_to_exteriorderivative(partial_derivatives(sfun, x))

# ── koszul ────────────────────────────────────────────────────────────────────
# κ_x(ω) = ι_x(ω) contracts with the evaluation point

κω1 = koszul(ω1)
@test evaluate(κω1, x) === TensorValues.koszul(x, ω1x)
@test evaluate(κω1, x) isa ExteriorFormValue{0,cov2}
# ι_{(1,1)}(3 dx¹ + 9 dx²) = 3 + 9
@test evaluate(κω1, x).data[1] ≈ 12.0   atol=1e-12

κvals = evaluate(κω1, xs)
@test length(κvals) == 3
# At (2,1): ω = (6,36), ι = 6*2 + 36*1
@test κvals[2].data[1] ≈ 48.0   atol=1e-12

ω_id = GenericField(y -> CoVectorValue(y[1], y[2]))
x0 = Point(2.0, 3.0)
# κ(ω)(x₀) = (x¹)² + (x²)² = 4 + 9
@test evaluate(koszul(ω_id), x0).data[1] ≈ 13.0   atol=1e-12
#@test evaluate(koszul(y -> ExteriorFormValue{1,Cov(2)}((y[1], y[2]))), x0).data[1] ≈ 13.0
@test evaluate(koszul(y -> CoVectorValue(y[1], y[2])), x0).data[1] ≈ 13.0

# κ of a 2-form: ι_{(2,3)}(dx¹∧dx²) = -x² dx¹ + x¹ dx²
ω_vol = GenericField(_ -> ExteriorFormValue{2,Cov(2)}((1.0,)))
κω_vol = evaluate(koszul(ω_vol), x0)
@test κω_vol isa CoVectorValue{cov2}
@test collect(κω_vol.data) ≈ [-3.0, 2.0]   atol=1e-12

# ── codifferential ────────────────────────────────────────────────────────────
# ω = 3x¹x² dx¹ + 5(x²)² dx²; ⋆ω = -5(x²)² dx¹ + 3x¹x² dx²
# d⋆ω = (∂₁(3x¹x²) + ∂₂(5(x²)²)) vol = (3x² + 10x²) vol, δω = -13x²

ωfun(x) = ExteriorFormValue{1,Cov(2)}((3.0*x[1]*x[2], 5.0*x[2]^2))
ωfld = GenericField(ωfun)
x1 = Point(1.0, 2.0)

δω = codifferential(ωfld,LinearAlgebra.I)
@test evaluate(δω, x1) === partials_to_flat_codifferential(partial_derivatives(ωfun, x1))
@test evaluate(δω, x1) isa ExteriorFormValue{0,cov2}
@test evaluate(δω, x1).data[1] ≈ -26.0   atol=1e-8

# ── pointwise Hodge star, and differentiating through it ──────────────────────
# The star is an `Operation`, and `partial_derivatives` differentiates the
# starred field as a black box, so the two compose.

star_ω = Operation(flat_hodge_star)(ωfld)
@test evaluate(star_ω, x1) === flat_hodge_star(evaluate(ωfld, x1))

d_star_ω = exterior_derivative(star_ω)
@test evaluate(d_star_ω, x1) ===
      partials_to_exteriorderivative(partial_derivatives(y -> flat_hodge_star(ωfun(y)), x1))
# δω = (-1)^{D(K-1)+1} ⋆d⋆ω, with D=2 and K=1 the sign is -1
@test (-flat_hodge_star(evaluate(d_star_ω, x1))).data[1] ≈ evaluate(δω, x1).data[1]   atol=1e-8

# ── lie_derivative ────────────────────────────────────────────────────────────

ufun(x) = ContraVectorValue(x[1] + 2.0*x[2], x[2]^2)
ufld = GenericField(ufun)

Lω = lie_derivative(ufld, ωfld)
@test evaluate(Lω, x1) === partials_to_liederivative(
  ufun(x1), ωfun(x1),
  partial_derivatives(ufun, x1), partial_derivatives(ωfun, x1))
@test evaluate(𝓛(ufld, ωfld), x1) === evaluate(Lω, x1)

# Broadcasting a two argument operator passes both of its arguments
bL = Broadcasting(lie_derivative)(Fill(ufld,3), Fill(ωfld,3))
@test evaluate(bL[1], x1) === evaluate(Lω, x1)

bd = Broadcasting(exterior_derivative)(Fill(ωfld,3))
@test evaluate(bd[1], x1) === evaluate(exterior_derivative(ωfld), x1)

bδ = Broadcasting(codifferential)(Fill(ωfld,3),LinearAlgebra.I)
@test evaluate(bδ[1], x1) === evaluate(δω, x1)

bκ = Broadcasting(koszul)(Fill(ωfld,3))
@test evaluate(bκ[1], x1) === evaluate(koszul(ωfld), x1)

end # module
