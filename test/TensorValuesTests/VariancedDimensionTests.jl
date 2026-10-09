module VariancedDimensionTests

using Test
using Gridap
using Gridap.TensorValues
using LinearAlgebra
using StaticArrays

t = TensorValue{con2, cov2}(1, 2, 3, 4)
@test 2 * t == TensorValue{Con(2), Cov(2)}(2, 4, 6, 8)
@test t + t == 2 * t

t2 = TensorValue{2, 2}(1:4...)
@test_throws "failed to change any argument" t + t2

v = VectorValue(1,1)
t = TensorValue{con2, cov2}(1, 1, 1, 1)
v_outer_t = ThirdOrderTensorValue{2, con2, cov2}(1,1,1,1,1,1,1,1)
@test v ⊗ t == v_outer_t
@test t ⊗ v != v_outer_t

ω13 = CoVectorValue{3}((1,2,3))
@test size(ω13) == (cov3,)
@test ω13 isa CoVectorValue{cov3}
@test ω13 ⊗ ω13 == TensorValue{cov3, cov3}(1, 2, 3, 2, 4, 6, 3, 6, 9)

ω23 = ExteriorFormValue{2,3}((1,2,3))
@test size(ω23) == (cov3,cov3)

# A fully symmetric tensor carries one variance per index, all of them being
# free as long as the dimensions are equal in value.
s3 = SymNTensorValue{Tuple{con3,con3,con3}}(1:10...)
@test size(s3) == (con3,con3,con3)
@test tensor_type(s3) == (3,0,0)
@test num_indep_components(s3) == 10

s_mix = SymNTensorValue{Tuple{cov2,con2}}(1,2,3)
@test size(s_mix) == (cov2,con2)
@test tensor_type(s_mix) == (1,1,0)
@test s_mix[1,2] == s_mix[2,1]

# Sizes equal in value but not in variance are distinct shapes: they neither
# compare equal nor promote together.
s_int = SymNTensorValue{Tuple{2,2}}(1,2,3)
s_cov = SymNTensorValue{Tuple{cov2,cov2}}(1,2,3)
@test s_int != s_cov
@test s_cov != s_mix
@test s_cov == SymNTensorValue{Tuple{cov2,cov2}}(1,2,3)
@test_throws "failed to change any argument" s_int + s_cov
@test_throws "failed to change any argument" s_cov + s_mix
@test promote_type(SymNTensorValue{Tuple{cov2,cov2},Int,2,3},
                   SymNTensorValue{Tuple{cov2,cov2},Float64,2,3}) ==
      SymNTensorValue{Tuple{cov2,cov2},Float64,2}

# Converting to an array drops the variance, as it does for any tensor
@test get_array(s_mix) isa SArray{Tuple{2,2},Int}
@test get_array(s_mix) == get_array(s_int)
@test Mutable(s_mix) == MArray{Tuple{2,2},Int,2,4}

@test tensor_type(1) == (0,0,0)
@test tensor_type(v) == (0,0,1)
@test tensor_type(t) == (1,1,0)
@test tensor_type(ω13) == (0,1,0)
@test tensor_type(ω23) == (0,2,0)
@test tensor_type(v_outer_t) == (1,1,1)

# ==================================================================
# Dimension predicates
# ==================================================================

@test is_same_variance(2, 2)
@test is_same_variance(cov2, cov3)
@test !is_same_variance(cov2, con2)
@test !is_same_variance(cov2, 2)

@test can_contract(2, 2)
@test !can_contract(2, 3)
@test can_contract(cov2, con2)
@test can_contract(con2, cov2)
@test !can_contract(cov2, cov2)
@test !can_contract(con2, con2)
@test !can_contract(cov2, 2)
@test !can_contract(cov2, con3)

@test inverse_variance(2) === 2
@test inverse_variance(cov2) === con2
@test inverse_variance(con2) === cov2
for D in (2, cov2, con2)
  @test can_contract(D, inverse_variance(D))
  @test inverse_variance(inverse_variance(D)) === D
end

# Negating a dimension preserves its variance
@test -cov2 === Cov(-2)
@test -con2 === Con(-2)

# ==================================================================
# Values used throughout
# ==================================================================

t_int   = TensorValue{2,2}(1,2,3,4)         # un-varianced, [1 3; 2 4]
t_mix   = TensorValue{con2,cov2}(1,2,3,4)   # (1,1) tensor
t_mix2  = TensorValue{cov2,con2}(1,2,3,4)   # (1,1) tensor, other way round
t_cov   = TensorValue{cov2,cov2}(1,2,3,4)   # (0,2) tensor
t_con   = TensorValue{con2,con2}(1,2,3,4)   # (2,0) tensor

s_int   = SymTensorValue{2}(1,2,3)          # [1 2; 2 3]
s_cov   = SymTensorValue{cov2}(1,2,3)
s_con   = SymTensorValue{con2}(1,2,3)

k_int   = SkewSymTensorValue{2}(1)          # [0 1; -1 0]
k_cov   = SkewSymTensorValue{cov2}(1)
k_con   = SkewSymTensorValue{con2}(1)

q_int   = SymTracelessTensorValue{2}(1,2)   # [1 2; 2 -1]
q_cov   = SymTracelessTensorValue{cov2}(1,2)
q_con   = SymTracelessTensorValue{con2}(1,2)

f_int   = SymFourthOrderTensorValue{2}(1:9...)
f_cov   = SymFourthOrderTensorValue{cov2}(1:9...)
f_con   = SymFourthOrderTensorValue{con2}(1:9...)

ω1      = CoVectorValue{cov2}((1,2))

# ==================================================================
# Comparison
# ==================================================================

@test t_mix ≈ t_mix
@test t_int ≈ TensorValue{2,2}(1.0,2.0,3.0,4.0)
# variance is part of the size, so tensors differing by it are not comparable
@test_throws DimensionMismatch (t_mix ≈ t_int)
@test_throws DimensionMismatch (t_mix ≈ t_cov)
@test_throws DimensionMismatch (t_mix ≈ t_mix2)
@test_throws DimensionMismatch (s_cov ≈ s_int)
# a genuine size mismatch is still reported by the array comparison
@test_throws DimensionMismatch (t_int ≈ TensorValue{2,3}(1:6...))

# Point approx with ContraVectorValue
@test VectorValue(1.0,2.0) ≈ ContraVectorValue{con2}((1.0,2.0))
@test ContraVectorValue{con2}((1.0,2.0)) ≈ VectorValue(1.0,2.0)
@test !(VectorValue(1.0,2.0) ≈ ContraVectorValue{con2}((1.0,3.0)))
@test_throws DimensionMismatch (VectorValue(1.0,2.0) ≈ CoVectorValue{cov2}((1.0,2.0)))
@test_throws DimensionMismatch (VectorValue(1.0,2.0) ≈ ContraVectorValue{con3}((1.0,2.0,3.0)))

# `==` already discriminates variance through the size type parameter
@test t_mix != t_int
@test t_mix != t_cov
@test t_mix == TensorValue{con2,cov2}(1.0,2.0,3.0,4.0)

# Comparison against a scalar covers every component, including those the `data`
# field does not store
@test !(0 <= k_int)
@test !isless(-1, k_int)
@test isless(-2, k_int)
@test 0 <= SkewSymTensorValue{2}(0)
@test !(0 <= q_int)
@test isless(-2, q_int)
# unchanged on tensors that store all their components
@test 0 <= VectorValue(1,3,3)
@test 1 <= VectorValue(1,3,3)
@test !isless(1, VectorValue(1,3,3))
@test isless(0, VectorValue(1,3,3))

# ==================================================================
# Promotion
# ==================================================================

# No promotion between tensors whose sizes differ in variance
@test_throws "failed to change any argument" t_mix + t_int
@test_throws "failed to change any argument" t_cov + t_con
@test_throws "failed to change any argument" s_cov + s_int
@test_throws "failed to change any argument" k_cov + k_int
@test_throws "failed to change any argument" f_cov + f_int
@test_throws "failed to change any argument" q_cov + q_con
@test_throws "failed to change any argument" ω1 + VectorValue(1,2)

# Promotion between the square tensor types survives a shared variance
@test t_cov + s_cov == TensorValue{cov2,cov2}(2,4,5,7)
@test (t_cov + s_cov) isa TensorValue{cov2,cov2}
@test t_cov + k_cov == TensorValue{cov2,cov2}(1,1,4,4)
@test (s_cov + q_cov) isa SymTensorValue{cov2}
# but not a mixed one: a (1,1) tensor is not of the type of the (0,2) tensors
@test_throws "failed to change any argument" t_mix + s_cov

# Same variance and same value is `===` on the size tuple
@test size(t_cov) === size(s_cov)
@test size(t_cov) !== size(t_con)
@test size(t_cov) == size(t_con)

# ==================================================================
# Conversion to plain arrays drops the variance
# ==================================================================

for a in (t_mix, t_cov, s_cov, k_cov, q_cov, ω1, ω23)
  @test size(get_array(a)) === Int.(size(a))
  @test size(SArray(a)) === size(get_array(a))
  @test size(mutable(a)) === size(get_array(a))
  @test all(i -> i isa Int, size(Mutable(typeof(a))))
end

@test get_array(t_mix) isa SMatrix{2,2,Int}
@test get_array(t_mix) == [1 3; 2 4]
@test Matrix(t_mix) == [1 3; 2 4]
@test Mutable(typeof(t_mix)) === MMatrix{2,2,Int,4}
@test SArray(ω23) isa SArray{Tuple{3,3},Int}
@test length(t_mix) === 4
@test num_components(typeof(t_mix)) === 4

# ==================================================================
# Matrix division
# ==================================================================

# Un-varianced behaviour is unchanged
@test one(TensorValue{3,3,Int}) \ VectorValue(1,2,3) == VectorValue(1,2,3)
@test one(SymTensorValue{3,Int}) \ VectorValue(1,2,3) == VectorValue(1,2,3)

# `inv(a)⋅b` is not symmetric even when `b` is, so the result must not be typed
# after `b`: a `SymTensorValue` result would silently drop the asymmetry.
a_ns = TensorValue{2,2}(1,2,3,5)
b_sy = SymTensorValue{2}(1,2,4)
@test !(a_ns \ b_sy isa AbstractSymTensorValue)
@test get_array(a_ns \ b_sy) ≈ get_array(a_ns) \ get_array(b_sy)
@test get_array(a_ns \ SkewSymTensorValue{2}(1)) ≈ get_array(a_ns) \ get_array(SkewSymTensorValue{2}(1))
# `a⋅x = b` round trip, which a symmetrised result would fail
@test a_ns ⋅ (a_ns \ b_sy) ≈ TensorValue{2,2}(get_array(b_sy))

# `a \ b` contracts the first index of `inv(a)` -- the second index of `a` with
# flipped variance -- against the first index of `b`, so `a` and `b` must share
# their first index.
g = TensorValue{con2,cov2}(2,0,0,4)
r = g \ TensorValue{con2,cov2}(2,0,0,4)
@test r isa TensorValue{con2,cov2}
@test r == TensorValue{con2,cov2}(1,0,0,1)

# first indices differ in variance
@test_throws DimensionMismatch g \ TensorValue{cov2,cov2}(1,2,3,4)
@test_throws DimensionMismatch g \ CoVectorValue{cov2}((1,2))
# left operand not square
@test_throws DimensionMismatch TensorValue{2,3}(1:6...) \ VectorValue(1,2)
# an order-0 right operand is not supported
@test_throws ErrorException t_int \ ExteriorFormValue{0,2}((1,))

# `inv` of a (0,2) tensor is a (2,0) tensor, so dividing a covariant vector by it
# gives a contravariant vector
q = TensorValue{cov2,cov2}(2,0,0,4) \ CoVectorValue{cov2}((1,2))
@test q isa ContraVectorValue{con2}
@test TensorValue{cov2,cov2}(2,0,0,4) ⋅ q ≈ CoVectorValue{cov2}((1.0,2.0))
# with a second index to carry, the result type does exist
r = TensorValue{cov2,cov2}(2,0,0,4) \ TensorValue{cov2,cov2}(2,0,0,4)
@test r isa TensorValue{con2,cov2}
@test r == TensorValue{con2,cov2}(1,0,0,1)

# `a / b` contracts the second indices, and flips the variance of `b`'s first
A = TensorValue{2,2}(1,2,3,5)
B = TensorValue{2,2}(2,0,1,3)
@test A / B ≈ A ⋅ inv(B)
@test (A / B) ⋅ B ≈ A
@test A / B isa TensorValue{2,2}

Acov = TensorValue{cov2,cov2}(1,2,3,5)
Bcov = TensorValue{cov2,cov2}(2,0,1,3)
@test (Acov / Bcov) isa TensorValue{cov2,con2}
@test get_array(Acov / Bcov) ≈ get_array(A / B)
# `x⋅b = a` round trip: the result contracts back against `b`
@test (Acov / Bcov) ⋅ Bcov ≈ Acov

@test_throws DimensionMismatch Acov / TensorValue{con2,con2}(2,0,1,3)
@test_throws DimensionMismatch TensorValue{2,3}(1:6...) / TensorValue{2,2}(1:4...)

# ==================================================================
# Inner product
# ==================================================================

# `inner` pairs index k of one operand with index k of the other, so the two
# must live in the same tensor space: equal sizes, variances included. A pair of
# opposite variances is a `⋅`, not a `⊙`.
@test_throws ErrorException s_cov ⊙ s_con
@test_throws ErrorException s_con ⊙ s_cov
@test_throws ErrorException q_cov ⊙ q_con
@test_throws ErrorException k_cov ⊙ k_con
@test_throws ErrorException f_cov ⊙ f_con
@test_throws ErrorException t_mix ⊙ t_mix2
@test_throws ErrorException s_cov ⊙ SymTensorValue{con3}(1:6...)
# operands of different order are refused before any index is paired
@test_throws ErrorException ω1 ⊙ t_cov
# same order, same variance, different dimensions
@test_throws ErrorException VectorValue(1,2,3) ⊙ VectorValue(1,2)
# `⋅` is what contracts a vector against a covector
@test contracted_product(Val(2), t_mix, t_mix2) == t_int ⊙ t_int

# Un-varianced indices contract on their own
@test s_int ⊙ s_int == 1 + 2*(2*2) + 3*3
@test t_int ⊙ t_int == 1 + 4 + 9 + 16

# A varianced index pairs with one of the same variance, which only the metric
# turns into an invariant scalar, so a metric-free call is refused
@test_throws ArgumentError s_cov ⊙ s_cov
@test_throws ArgumentError s_con ⊙ s_con
@test_throws ArgumentError q_cov ⊙ q_cov
@test_throws ArgumentError k_cov ⊙ k_cov
@test_throws ArgumentError f_cov ⊙ f_cov
@test_throws ArgumentError t_cov ⊙ t_cov
@test_throws ArgumentError t_mix ⊙ t_mix
@test_throws ArgumentError ω1 ⊙ ω1

# The metric is passed as to every metric dependent operation: either tensor,
# in any order, `LinearAlgebra.I` for the Euclidean product
gm_lo = SymTensorValue{cov2}(2.0, 0.5, 3.0)
gm_up = inv(gm_lo)
ω1f   = CoVectorValue{cov2}((1.0, 2.0))

@test inner(s_cov, s_cov, gm_lo) ≈ inner(s_cov, s_cov, gm_up)
@test inner(s_cov, s_cov, I) == flat_inner(s_cov, s_cov) == s_int ⊙ s_int
@test inner(t_mix, t_mix, I) == t_int ⊙ t_int
@test_throws DimensionMismatch inner(s_cov, s_cov, SymTensorValue{cov3}(1:6...))

# A covariant pair contracts through `g_inv` and a contravariant one through
# `g`, as `alt_inner` does, which is what makes the scalar invariant: raising
# both operands does not change it
@test inner(ω1f, ω1f, gm_lo) ≈ sum(ω1f[i]*gm_up[i,j]*ω1f[j] for i in 1:2, j in 1:2)
@test inner(ω1f, ω1f, gm_lo) ≈ inner(sharp(ω1f, gm_up), sharp(ω1f, gm_up), gm_lo)

# One metric cannot serve indices of different dimensions
@test_throws DimensionMismatch inner(TensorValue{cov2,cov3}(1.:6...),
                                     TensorValue{cov2,cov3}(1.:6...), gm_lo)

# Every operation built on `inner` takes the same metric arguments and forwards
# them, which is what generalizes it to varianced tensors
@test norm(ω1f, gm_lo) ≈ sqrt(inner(ω1f, ω1f, gm_lo))
@test meas(ω1f, gm_lo) == norm(ω1f, gm_lo)
@test norm(t_cov, gm_lo) ≈ sqrt(inner(t_cov, t_cov, gm_lo))
@test norm(ω1f, I) == norm(VectorValue(1.0, 2.0))

# the exponent of `norm(u, p)` comes before the metric arguments, and only the
# 2-norm of a tensor is defined
@test norm(VectorValue(3.0, 4.0), 2) == 5.0
@test norm(ω1f, 2, gm_lo) == norm(ω1f, gm_lo)
@test norm(t_cov, 2, gm_lo) == norm(t_cov, gm_lo)
@test_throws ArgumentError norm(VectorValue(3.0, 4.0), 1)
@test_throws ArgumentError norm(ω1f, Inf, gm_lo)

# the dual of a component basis is read off the inverse of its Gram matrix, so
# it needs the metric of a varianced type
Vcov = SymTensorValue{cov2,Float64}
@test_throws ArgumentError representatives_of_componentbasis_dual(Vcov)
basis_cov = component_basis(Vcov)
dual_cov  = representatives_of_componentbasis_dual(Vcov, gm_lo)
@test [ inner(bi, dj, gm_lo) for bi in basis_cov, dj in dual_cov ] ≈ I
@test representatives_of_basis_dual(basis_cov, gm_lo) ≈ dual_cov

# A skew symmetric tensor is orthogonal to every symmetric one under any
# symmetric metric, so that pair needs none
@test iszero(k_int ⊙ s_int)
@test iszero(k_cov ⊙ s_cov)
@test iszero(s_cov ⊙ k_cov)
@test inner(k_cov, s_cov, gm_lo) ≈ 0 atol=1e-14
@test_throws ErrorException k_cov ⊙ s_con

# `inner` of an order-4 with an order-2 tensor is `double_contraction`, through
# a deprecated pair of methods that pair the indices by rank and not by
# position. They inherit the variance check of `double_contraction`, so two
# covariant arguments do not contract
@test get_array(f_int ⊙ s_int) == get_array(contracted_product(Val(2), f_int, s_int))
@test get_array(s_int ⊙ f_int) == get_array(contracted_product(Val(2), s_int, f_int))
@test_throws DimensionMismatch f_cov ⊙ s_cov
@test_throws DimensionMismatch s_cov ⊙ f_cov
@test_throws ErrorException f_cov ⊙ s_con
@test_throws ErrorException f_int ⊙ TensorValue{2,3}(1:6...)

# ==================================================================
# Double contraction
# ==================================================================

@test get_array(f_cov ⋅² s_con) == get_array(f_int ⋅² s_int)
@test (f_cov ⋅² s_con) isa SymTensorValue{cov2}
@test get_array(s_cov ⋅² f_con) == get_array(s_int ⋅² f_int)
@test (s_cov ⋅² f_con) isa SymTensorValue{con2}
@test_throws DimensionMismatch f_cov ⋅² s_cov
@test_throws DimensionMismatch s_cov ⋅² f_cov

# Non symmetric second order operands are symmetrised first, which requires
# their two indices to have the same variance
@test get_array(f_cov ⋅² t_con) == get_array(f_int ⋅² t_int)
@test get_array(t_cov ⋅² f_con) == get_array(t_int ⋅² f_int)
@test_throws ErrorException f_cov ⋅² t_mix
@test_throws ErrorException t_mix ⋅² f_cov

# Third order against fourth order symmetric
r3  = ThirdOrderTensorValue{2,cov2,cov2}(1:8...)
r3i = ThirdOrderTensorValue{2,2,2}(1:8...)
@test get_array(r3 ⋅² f_con) == get_array(r3i ⋅² f_int)
@test (r3 ⋅² f_con) isa ThirdOrderTensorValue{2,con2,con2}
@test_throws DimensionMismatch r3 ⋅² f_cov

# Two fourth order symmetric tensors of opposite variance contract to a tensor
# with two indices of each variance, which no `SymFourthOrderTensorValue` holds
@test_throws DimensionMismatch f_cov ⋅² f_cov
@test (f_int ⋅² f_int) isa SymFourthOrderTensorValue{2}
@test (f_cov ⋅² f_con) isa NTensorValue{Tuple{cov2,cov2,con2,con2}}
@test get_array(f_cov ⋅² f_con) == get_array(f_int ⋅² f_int)

# ==================================================================
# Congruent product
# ==================================================================

# bᵀ⋅a⋅b : the first index of `b` is contracted against the indices of `a`
b_con = TensorValue{con2,con2}(1,2,3,4)
b_int = TensorValue{2,2}(1,2,3,4)
@test get_array(congruent_prod(t_cov, b_con)) == get_array(congruent_prod(t_int, b_int))
@test congruent_prod(t_cov, b_con) isa TensorValue{con2,con2}
@test congruent_prod(s_cov, b_con) isa SymTensorValue{con2}
@test congruent_prod(k_cov, b_con) isa SkewSymTensorValue{con2}
@test_throws DimensionMismatch congruent_prod(t_cov, TensorValue{cov2,cov2}(1,2,3,4))

# ==================================================================
# Determinant, inverse, eigen, sqrt
# ==================================================================

@test det(t_mix)  == det(t_int)
@test det(t_mix2) == det(t_int)
@test det(t_cov) == det(t_int)
@test det(t_con) == det(t_int)
@test det(s_cov) == det(s_int)
@test det(k_cov) == det(k_int)
@test det(SkewSymTensorValue{3}(1,2,3)) == 0
@test det(SkewSymTensorValue{cov3}(1,2,3)) == 0

e_mix = TensorValue{con2,cov2}(2.,1.,1.,3.)
e_int = TensorValue{2,2}(2.,1.,1.,3.)
@test eigen(e_mix).values ≈ eigen(e_int).values
@test_throws ErrorException eigen(TensorValue{cov2,cov2}(2.,1.,1.,3.))

# `inv` raises the lower indices and lowers the upper ones
@test inv(t_int) isa TensorValue{2,2}
@test inv(t_mix) isa TensorValue{con2,cov2}
@test inv(t_mix2) isa TensorValue{cov2,con2}
@test inv(t_cov) isa TensorValue{con2,con2}
@test inv(t_con) isa TensorValue{cov2,cov2}
@test get_array(inv(t_cov)) ≈ get_array(inv(t_int))
@test inv(s_cov) isa SymTensorValue{con2}
@test inv(k_cov) isa SkewSymTensorValue{con2}
@test get_array(inv(s_cov)) ≈ get_array(inv(s_int))
@test_throws ErrorException inv(TensorValue{2,3}(1:6...))
# a (1,1) tensor composes with its inverse into the identity of the same shape
@test inv(t_mix) ⋅ t_mix ≈ TensorValue{con2,cov2}(1.,0.,0.,1.)

# `sqrt` keeps the indices of its argument
@test sqrt(TensorValue{2,2}(4.,0.,0.,4.)) isa TensorValue{2,2}
@test sqrt(TensorValue{con2,cov2}(4.,0.,0.,4.)) isa TensorValue{con2,cov2}
@test sqrt(TensorValue{con2,cov2}(4.,0.,0.,4.)) == TensorValue{con2,cov2}(2.,0.,0.,2.)
@test_throws ErrorException sqrt(TensorValue{cov2,cov2}(4.,0.,0.,4.))

# ==================================================================
# Trace
# ==================================================================

@test tr(t_mix) == tr(t_int)
@test tr(t_mix2) == tr(t_int)
@test_throws ArgumentError tr(t_cov)
@test_throws ArgumentError tr(t_con)
@test_throws ArgumentError tr(s_cov)
@test_throws ArgumentError tr(q_cov)
@test_throws ArgumentError tr(k_cov)
@test iszero(tr(q_int))
@test iszero(tr(k_int))
@test_throws ArgumentError tr(TensorValue{2,3}(1:6...))

# Third order: only the two traced indices are constrained
h3  = ThirdOrderTensorValue{con2,cov2,cov3}(1:12...)
h3i = ThirdOrderTensorValue{2,2,3}(1:12...)
@test tr(h3i) isa VectorValue{3}
@test tr(h3) isa CoVectorValue{cov3}
@test get_array(tr(h3)) == get_array(tr(h3i))
@test_throws ArgumentError tr(ThirdOrderTensorValue{cov2,cov2,3}(1:12...))

# `tr` and the self-contraction agree, variance included
@test tr(h3) == tensor_contraction(h3, Val(1), Val(2))
@test tr(t_mix) == tensor_contraction(t_mix, Val(1), Val(2))
@test_throws DimensionMismatch tensor_contraction(t_cov, Val(1), Val(2))

# ==================================================================
# meas and norm
# ==================================================================

@test_throws ErrorException meas(t_mix)
@test_throws ErrorException meas(t_cov)
# a row read as a vector is only length preserving without variance
@test meas(TensorValue{1,2}(3,4)) == 5
@test_throws ErrorException meas(TensorValue{1,cov2}(3,4))
# a non-square Jacobian measures the span of its rows, a point domain having the
# unit measure that leaves an integral summing values
@test meas(TensorValue{2,3}(2,0, 0,3, 0,0)) == 6
@test meas(TensorValue{3,4}(1,0,0, 0,1,0, 0,0,1, 0,0,0)) == 1
@test meas(zero(TensorValue{0,2,Float64,0})) == 1

# the transposed Jacobian of a chart scales that one chart's volume form, so it
# is measured with its variance dropped, at any shape and only in that order
@test meas(TensorValue{cov2,con2}(1,2,3,4)) == meas(t_int)
@test meas(TensorValue{cov2,con3}(2,0, 0,3, 0,0)) == 6
@test meas(zero(TensorValue{Cov(0),con2,Float64,0})) == 1
@test_throws ErrorException meas(TensorValue{con2,cov3}(1,0, 0,1, 0,0))

# `norm` inherits its variance rules from `inner`, which pairs index k of `u`
# with index k of `u`: no varianced tensor has an invariant Frobenius norm, a
# (1,1) tensor included -- the invariant quantity there is `tr(u⋅u)`.
@test norm(t_int) == sqrt(1 + 4 + 9 + 16)
@test_throws ArgumentError norm(t_cov)
@test_throws ArgumentError norm(t_mix)
@test tr(t_mix ⋅ t_mix) == tr(t_int ⋅ t_int)
# a form carries covariant indices, so its Euclidean norm goes through
# `flat_inner`, which raises them by the identity metric
@test_throws ArgumentError norm(ω1)
@test sqrt(flat_inner(ω1, ω1)) == norm(VectorValue(1,2))
@test ω1 ⊙ᴱ ω1 === flat_inner(ω1, ω1)

# ==================================================================
# Variance of a metric tensor
# ==================================================================

# A metric lowers indices and is a (0,2) tensor; its inverse raises them and is
# a (2,0) one. Its variance is what says which of the two it is, so an
# un-varianced one is rejected.
ω2    = CoVectorValue{cov2}((1.0, 2.0))
v2    = ContraVectorValue{con2}((1.0, 2.0))
v2_un = VectorValue(1.0, 2.0)
g_id  = SymTensorValue{2}(1.0, 0.0, 1.0)
g_lo  = SymTensorValue{cov2}(1.0, 0.0, 1.0)
g_up  = SymTensorValue{con2}(1.0, 0.0, 1.0)

@test tensor_type(g_lo) == (0,2,0)
@test tensor_type(g_up) == (2,0,0)
@test tensor_type(g_id) == (0,0,2)

# `hodge_star`, `sharp` and `alt_inner` raise indices with `g_inv`; `flat`
# lowers one with `g`
@test hodge_star(ω2, g_up, 1.0) ≈ flat_hodge_star(ω2)
@test sharp(ω2, g_up)           ≈ v2
@test alt_inner(ω2, ω2, g_up)   ≈ flat_alt_inner(ω2, ω2)
@test flat(v2, g_lo)            ≈ ω2

# Either tensor may be given, the other one is inverted on demand
@test hodge_star(ω2, g_lo, 1.0) ≈ hodge_star(ω2, g_up, 1.0)
@test sharp(ω2, g_lo)           ≈ sharp(ω2, g_up)
@test alt_inner(ω2, ω2, g_lo)   ≈ alt_inner(ω2, ω2, g_up)
@test flat(v2, g_up)            ≈ flat(v2, g_lo)

# √|det(g)| is computed from whichever tensor is given when it is omitted
@test hodge_star(ω2, g_up) ≈ hodge_star(ω2, g_up, 1.0)
@test hodge_star(ω2, g_lo) ≈ hodge_star(ω2, g_up, 1.0)

# An un-varianced metric does not say whether it raises or lowers
@test_throws ArgumentError hodge_star(ω2, g_id, 1.0)
@test_throws ArgumentError sharp(ω2, g_id)
@test_throws ArgumentError alt_inner(ω2, ω2, g_id)
@test_throws ArgumentError flat(v2, g_id)

# `alt_inner` also accepts a K-vector pair, contracting through `g` instead of
# `g_inv`, symmetrically to how `flat`/`sharp` lower/raise a vector/form
@test alt_inner(v2, v2, g_lo) ≈ flat_alt_inner(v2, v2)
@test alt_inner(v2, v2, g_up) ≈ alt_inner(v2, v2, g_lo)
@test_throws ArgumentError alt_inner(v2, v2, g_id)

# A `VectorValue` is taken as contravariant, but the result still carries a variance
@test flat(v2_un, g_lo) ≈ flat(v2, g_lo)
@test flat(v2_un, g_lo) isa CoVectorValue{cov2}

# metric errors
@test_throws ArgumentError sharp(v2, g_up)
@test_throws ArgumentError flat(ω2, g_lo)
@test_throws DimensionMismatch sharp(ω2, SymTensorValue{con3}(1.,0.,0.,1.,0.,1.))
@test_throws DimensionMismatch flat(v2, SymTensorValue{cov3}(1.,0.,0.,1.,0.,1.))

# `LinearAlgebra.I` (a `UniformScaling`) selects the Euclidean (identity
# metric) version of each operator, without needing a `SymTensorValue`
@test hodge_star(ω2, I)     ≈ flat_hodge_star(ω2)
@test flat(v2, I)           ≈ flat(v2, g_lo)
@test sharp(ω2, I)          ≈ sharp(ω2, g_up)
@test alt_inner(ω2, ω2, I)  ≈ alt_inner(ω2, ω2, g_up)
@test alt_inner(v2, v2, I)  ≈ alt_inner(v2, v2, g_lo)
@test flat(v2, I) isa CoVectorValue{cov2}
@test sharp(ω2, I) isa ContraVectorValue{con2}

# `flat_flat`/`flat_sharp` are the Euclidean musical isomorphisms directly,
# inverses of each other, and reject the wrong variance like `flat`/`sharp` do
@test flat_flat(v2)  ≈ ω2
@test flat_sharp(ω2) ≈ v2
@test flat_sharp(flat_flat(v2)) === v2
@test flat_flat(flat_sharp(ω2)) === ω2
@test flat_flat(v2_un) ≈ flat_flat(v2)
@test_throws ArgumentError flat_flat(ω2)
@test_throws ArgumentError flat_sharp(v2)

# the flat operators alias the non-flat symbol suffixed by ᴱ, like flat_hodge_star/⋆ᴱ
@test ♭ᴱ(v2) === flat_flat(v2)
@test ♯ᴱ(ω2) === flat_sharp(ω2)

# `flat`/`sharp` lower/raise one chosen index of a tensor of any order, the
# index being given by `Val`; the other indices keep their variance
m_up = TensorValue{con2,con2}(1.0, 2.0, 3.0, 4.0)
m_lo = TensorValue{cov2,cov2}(1.0, 2.0, 3.0, 4.0)
@test flat(m_up, Val(1), g_lo) isa TensorValue{cov2,con2}
@test flat(m_up, Val(2), g_lo) isa TensorValue{con2,cov2}
@test sharp(m_lo, Val(1), g_up) isa TensorValue{con2,cov2}
@test flat(m_up, Val(1), gm_lo) ≈ TensorValue{cov2,con2}(
  (sum(gm_lo[i,l]*m_up[l,j] for l in 1:2) for i in 1:2, j in 1:2)...)
@test flat(m_up, Val(2), gm_lo) ≈ TensorValue{con2,cov2}(
  (sum(gm_lo[j,l]*m_up[i,l] for l in 1:2) for i in 1:2, j in 1:2)...)
@test sharp(m_lo, Val(2), gm_up) ≈ TensorValue{cov2,con2}(
  (sum(gm_up[j,l]*m_lo[i,l] for l in 1:2) for i in 1:2, j in 1:2)...)

# the two are inverse of each other, index by index
@test sharp(flat(m_up, Val(1), gm_lo), Val(1), gm_up) ≈ m_up
@test flat(sharp(m_lo, Val(2), gm_up), Val(2), gm_lo) ≈ m_lo

# lowering one index of a symmetric tensor gives a tensor that is not symmetric,
# so `flat` refuses to do it unless `drop_sym` allows it to drop the symmetry
@test_throws "does not preserve the symmetries" flat(
  SymTensorValue{con2}(1.0, 2.0, 3.0), Val(1), gm_lo)
@test flat(SymTensorValue{con2}(1.0, 2.0, 3.0), Val(1), gm_lo; drop_sym=Val(true)) isa
      TensorValue{cov2,con2}

# the indexed form agrees with the musical isomorphisms of a 1-tensor
@test flat(v2, Val(1), gm_lo)  ≈ flat(v2, gm_lo)
@test sharp(ω2, Val(1), gm_up) ≈ sharp(ω2, gm_up)

# `LinearAlgebra.I` selects the Euclidean version, which only relabels variance
@test flat(m_up, Val(1), I)  === flat_flat(m_up, Val(1))
@test sharp(m_lo, Val(2), I) === flat_sharp(m_lo, Val(2))
@test get_array(flat_flat(m_up, Val(1))) == get_array(m_up)

# an index can only be lowered if it is contravariant, and raised if covariant
@test_throws ArgumentError flat(m_lo, Val(1), g_lo)
@test_throws ArgumentError sharp(m_up, Val(1), g_up)
@test_throws ArgumentError flat(t_int, Val(1), g_lo)
@test_throws ArgumentError flat_flat(m_lo, Val(1))
@test_throws ArgumentError flat_sharp(m_up, Val(1))
# and only if it exists
@test_throws ArgumentError flat(m_up, Val(3), g_lo)
@test_throws ArgumentError flat(m_up, Val(0), g_lo)
@test_throws DimensionMismatch flat(m_up, Val(1), SymTensorValue{cov3}(1.,0.,0.,1.,0.,1.))

# several indices are lowered or raised in one call when the `Val` holds a list
# of them, which gives what doing them one after the other does
@test flat(m_up, Val((1,2)), gm_lo) ≈ flat(flat(m_up, Val(1), gm_lo), Val(2), gm_lo)
@test flat(m_up, Val((1,2)), gm_lo) isa TensorValue{cov2,cov2}
@test flat(m_up, Val((2,1)), gm_lo) ≈ flat(m_up, Val((1,2)), gm_lo)
@test flat(m_up, Val((1,)), gm_lo) === flat(m_up, Val(1), gm_lo)
@test sharp(m_lo, Val((1,2)), gm_up) ≈ sharp(sharp(m_lo, Val(1), gm_up), Val(2), gm_up)
@test flat(m_up, Val((1,2)), I) === flat_flat(m_up, Val((1,2)))

# an empty list of indices leaves the tensor alone
@test flat(m_up, Val(()), gm_lo) === m_up
@test sharp(m_lo, Val(()), gm_up) === m_lo

# each listed index must exist, be listed once, and have the right variance
@test_throws ArgumentError flat(m_up, Val((1,1)), gm_lo)
@test_throws ArgumentError flat(m_up, Val((1,3)), gm_lo)
@test_throws ArgumentError flat(t_mix, Val((1,2)), gm_lo)
# and one metric only acts on indices that share their dimension
@test_throws DimensionMismatch flat(TensorValue{con2,con3}(1.:6...), Val((1,2)), gm_lo)

# `inner` is those lowerings and raisings followed by a component-wise sum
# against the other operand
@test inner(m_up, m_up, gm_lo) ≈
  sum(get_array(m_up) .* get_array(flat(m_up, Val((1,2)), gm_lo)))

# A tensor type that encodes index symmetries keeps them when every index is
# flipped by the same metric, that is when the map is a congruence `G…G a G…G`.
# The result then holds only the independent components, so it is exactly
# symmetric and not symmetric up to roundoff.
s_up = SymTensorValue{con2}(1.0, 2.0, 3.0)
@test flat(s_up, Val((1,2)), gm_lo) isa SymTensorValue{cov2}
@test get_array(flat(s_up, Val((1,2)), gm_lo)) ≈
      get_array(gm_lo)*get_array(s_up)*get_array(gm_lo)'
let A = get_array(flat(s_up, Val((1,2)), gm_lo))
  @test A[1,2] === A[2,1]
end
# the order the indices are listed in does not matter, they are flipped at once
@test flat(s_up, Val((2,1)), gm_lo) === flat(s_up, Val((1,2)), gm_lo)
@test sharp(flat(s_up, Val((1,2)), gm_lo), Val((1,2)), gm_up) ≈ s_up

# a `SymTracelessTensorValue` gives a `SymTensorValue`: the congruence preserves
# the symmetry, but the Euclidean trace is not invariant under it
q_up = SymTracelessTensorValue{con2}(1.0, 2.0)
@test flat(q_up, Val((1,2)), gm_lo) isa SymTensorValue{cov2}
@test get_array(flat(q_up, Val((1,2)), gm_lo)) ≈
      get_array(gm_lo)*get_array(q_up)*get_array(gm_lo)'
# even Euclideanly, where the trace does survive: one policy for both metrics
@test flat_flat(q_up, Val((1,2))) isa SymTensorValue{cov2}

# the same holds for the fully symmetric tensors of any order
a3_up = SymNTensorValue{Tuple{con2,con2,con2}}(1.0, 2.0, 3.0, 4.0)
@test flat(a3_up, Val((1,2,3)), gm_lo) isa SymNTensorValue{Tuple{cov2,cov2,cov2}}
@test get_array(flat(a3_up, Val((1,2,3)), gm_lo)) ≈ [
  sum(gm_lo[i,p]*gm_lo[j,q]*gm_lo[k,r]*a3_up[p,q,r] for p in 1:2, q in 1:2, r in 1:2)
  for i in 1:2, j in 1:2, k in 1:2]
# and for the two minor symmetries of a `SymFourthOrderTensorValue`
c4_up = SymFourthOrderTensorValue{con2}(1.0:9.0...)
@test flat(c4_up, Val((1,2,3,4)), gm_lo) isa SymFourthOrderTensorValue{cov2}
@test get_array(flat(c4_up, Val((1,2,3,4)), gm_lo)) ≈ [
  sum(gm_lo[i,p]*gm_lo[j,q]*gm_lo[k,r]*gm_lo[m,s]*c4_up[p,q,r,s]
      for p in 1:2, q in 1:2, r in 1:2, s in 1:2)
  for i in 1:2, j in 1:2, k in 1:2, m in 1:2]

# flipping only part of the indices leaves a symmetry that no single tensor type
# carries, so it is refused rather than silently dropped
@test_throws "does not preserve the symmetries" flat(a3_up, Val((1,2)), gm_lo)
@test_throws "does not preserve the symmetries" flat(c4_up, Val((1,2)), gm_lo)
@test_throws "does not preserve the symmetries" flat_flat(s_up, Val(1))
@test_throws "does not preserve the symmetries" flat_sharp(
  SymTensorValue{cov2}(1.0, 2.0, 3.0), Val(1))

# the storage symmetry of a mixed variance tensor is not a tensorial property:
# flipping both indices of `a_i{}^j` is the similarity `g⋅A⋅g_inv`, which is
# symmetric in no pair of indices, so no index list preserves the symmetry
am = SymNTensorValue{Tuple{cov2,con2}}(1.0, 2.0, 3.0)
@test_throws "does not preserve the symmetries" flat(am, Val((2,)), gm_lo)
@test_throws "does not preserve the symmetries" sharp(am, Val((1,)), gm_up)

# `drop_sym=Val(true)` lifts the refusal and flips the indices of the full rank
# tensor, which is what a caller flipping the indices of an unknown tensor needs
@test flat(a3_up, Val((1,2)), gm_lo; drop_sym=Val(true)) isa
      ThirdOrderTensorValue{cov2,cov2,con2}
@test flat(am, Val((2,)), gm_lo; drop_sym=Val(true)) isa TensorValue{cov2,cov2}
@test flat_flat(s_up, Val(1); drop_sym=Val(true)) isa TensorValue{cov2,con2}
# an index list that does preserve the symmetries is honoured either way, so
# asking to drop them costs a caller nothing when there is nothing to drop
@test flat(s_up, Val((1,2)), gm_lo; drop_sym=Val(true)) ===
      flat(s_up, Val((1,2)), gm_lo)
@test flat(a3_up, Val((1,2,3)), gm_lo; drop_sym=Val(true)) ===
      flat(a3_up, Val((1,2,3)), gm_lo)
# a list the operation refuses for its variance still reports that, not symmetry
@test_throws "not all contravariant" flat(
  SymTensorValue{cov2}(1.0, 2.0, 3.0), Val((1,2)), gm_lo)

# an empty list leaves the tensor alone, symmetries included
@test flat(s_up, Val(()), gm_lo) === s_up
@test flat(s_up, Val(()), gm_lo; drop_sym=Val(true)) === s_up

# `inner` flips the indices of its second operand, which it knows nothing about,
# so it allows the symmetries to be dropped, and agrees with the dense sum
@test inner(s_up, s_up, gm_lo) ≈
  sum(get_array(s_up) .* get_array(flat(s_up, Val((1,2)), gm_lo)))
@test inner(q_up, q_up, gm_lo) ≈
  sum(get_array(q_up) .* get_array(flat(q_up, Val((1,2)), gm_lo)))
@test inner(am, am, gm_lo) ≈
  sum(get_array(am) .* get_array(sharp(flat(am, Val((2,)), gm_lo; drop_sym=Val(true)),
                                       Val((1,)), gm_up; drop_sym=Val(true))))

# `flat`/`sharp` of a skew-symmetric tensor act on all of its indices at once,
# each component of the result being summed over the sorted multi-indices with
# the determinant of the metric's minor on the two of them as weight
ω2f = ExteriorFormValue{2,cov3}((1.0, 2.0, 3.0))
v2f = ExteriorVectorValue{2,con3}((4.0, 5.0, 6.0))
g3  = SymTensorValue{cov3}(2.0, 0.3, 0.1, 2.5, 0.2, 3.0)
g3i = inv(g3)
@test sharp(ω2f, g3) isa ExteriorVectorValue{2,con3}
@test flat(v2f, g3)  isa ExteriorFormValue{2,cov3}
let c = ((1,2), (1,3), (2,3))
  minor(M, I, J) = M[I[1],J[1]]*M[I[2],J[2]] - M[I[1],J[2]]*M[I[2],J[1]]
  @test sharp(ω2f, g3) ≈ ExteriorVectorValue{2,con3}(
    ntuple(i -> sum(minor(g3i, c[i], J)*ω2f[J...] for J in c), 3))
  @test flat(v2f, g3) ≈ ExteriorFormValue{2,cov3}(
    ntuple(i -> sum(minor(g3, c[i], J)*v2f[J...] for J in c), 3))
end

# raising the single component of a top form scales it by det(g_inv)
@test sharp(ExteriorFormValue{3,cov3}((2.0,)), g3) ≈
      ExteriorVectorValue{3,con3}((2.0*det(get_array(g3i)),))

# the two are inverse of each other, and `I` gives the Euclidean version, which
# only relabels the variance
@test flat(sharp(ω2f, g3), g3) ≈ ω2f
@test sharp(ω2f, I) === flat_sharp(ω2f)
@test flat(v2f, I)  === flat_flat(v2f)
@test flat_sharp(ω2f) === ExteriorVectorValue{2,con3}(ω2f.data)

# lowering or raising only part of the indices of a skew-symmetric tensor would
# not leave it skew-symmetric, so the only index list it accepts is all of them
@test sharp(ω2f, Val((1,2)), g3) === sharp(ω2f, g3)
@test flat(v2f, Val((1,2)), g3)  === flat(v2f, g3)
@test flat_sharp(ω2f, Val((1,2))) === flat_sharp(ω2f)
@test_throws ArgumentError sharp(ω2f, Val((1,)), g3)
@test_throws ArgumentError sharp(ω2f, Val((2,1)), g3)
@test_throws ArgumentError flat(v2f, Val(1), g3)
@test_throws ArgumentError flat_sharp(ω2f, Val((1,)))
# unless `drop_sym` allows the full rank result, as for the symmetric types
@test sharp(ω2f, Val((1,)), g3; drop_sym=Val(true)) isa TensorValue{con3,cov3}
@test flat(v2f, Val(1), g3; drop_sym=Val(true))      isa TensorValue{cov3,con3}
@test flat_sharp(ω2f, Val((1,)); drop_sym=Val(true)) isa TensorValue{con3,cov3}

# a K-form contracts with itself through that one raising: the scalar agrees
# with the dense sum over the expanded tensors, and is the `K!` multiple of the
# alternating inner product that stays within the independent components
@test inner(ω2f, ω2f, g3) ≈ sum(get_array(ω2f) .* get_array(sharp(ω2f, g3)))
@test inner(v2f, v2f, g3) ≈ sum(get_array(v2f) .* get_array(flat(v2f, g3)))
@test inner(ω2f, ω2f, g3) ≈ factorial(2) * alt_inner(ω2f, ω2f, g3)
@test inner(v2f, v2f, g3) ≈ factorial(2) * alt_inner(v2f, v2f, g3)
@test inner(ω2f, ω2f, I) == flat_inner(ω2f, ω2f)
# a varianced pair refuses to guess a metric, as `alt_inner` itself does
@test_throws ArgumentError inner(ω2f, ω2f)
@test_throws ArgumentError alt_inner(ω2f, ω2f)

# `hodge_star` keeps the variance of `ω`: it raises the indices of a form with
# `g_inv` and scales by √|det(g)|, and does the dual for a K-vector
gd  = SymTensorValue{cov2,Float64,3}(2.0, 0.0, 3.0)
gdi = inv(gd)
sdg = sqrt(6.0)
for v in (ContraVectorValue{con2}((1.0, 0.0)), ContraVectorValue{con2}((0.7, -0.3)))
  star_v = hodge_star(v, gd)
  @test star_v isa ContraVectorValue{con2}
  # agrees with the form-side star transported by the musical isomorphisms
  @test star_v ≈ sharp(hodge_star(flat(v, gd), gdi), gdi)
  # ⋆⋆ = (-1)^{K(D-K)} = -1 for K = 1, D = 2
  @test hodge_star(star_v, gd) ≈ -1.0 * v
  # either tensor, and √|det(g)| given or not
  @test hodge_star(v, gdi)      ≈ star_v
  @test hodge_star(v, gd,  sdg) ≈ star_v
  @test hodge_star(v, gdi, sdg) ≈ star_v
end

# The volume factor of a K-vector is the reciprocal of a K-form's
@test hodge_star(ExteriorVectorValue{0,con2}((5.0,)), gd) ≈ ExteriorVectorValue{2,con2}((5.0/sdg,))
@test hodge_star(ExteriorFormValue{0,cov2}((5.0,)),   gd) ≈ ExteriorFormValue{2,cov2}((5.0*sdg,))

# An un-varianced K-tensor says neither, `flat_hodge_star` is its Euclidean star
@test_throws ErrorException hodge_star(SkewSymKTensorValue{Int,1,2}((1.0, 2.0)), gd)

# Each result carries the variance its own index ends up with
@test flat(v2, g_lo)            isa CoVectorValue{cov2}
@test hodge_star(ω2, g_up, 1.0) isa CoVectorValue{cov2}
@test sharp(ω2, g_up)           isa ContraVectorValue{con2}

# ==================================================================
# Adjoint, transpose and permutedims
# ==================================================================

# Reversing the slots carries each index along with its variance
@test transpose(t_int) == TensorValue{2,2}(1,3,2,4)
@test transpose(t_cov) isa TensorValue{cov2,cov2}
@test transpose(t_cov) == TensorValue{cov2,cov2}(1,3,2,4)
@test transpose(t_mix) isa TensorValue{cov2,con2}
@test transpose(t_mix) == TensorValue{cov2,con2}(1,3,2,4)
@test adjoint(t_mix) == transpose(t_mix)
@test transpose(transpose(t_mix)) === t_mix
@test transpose(TensorValue{cov2,cov3}(1:6...)) isa TensorValue{cov3,cov2}
@test transpose(TensorValue{cov2,3}(1:6...)) isa TensorValue{3,cov2}
# adjoint conjugates, which is orthogonal to the variance
@test adjoint(TensorValue{con2,cov2}(1+im,2,3,4)) == TensorValue{cov2,con2}(1-im,3,2,4)
# the symmetric types carry twice the same index, so they always transpose
@test transpose(s_cov) === s_cov
@test transpose(k_cov) == -k_cov

# Any permutation is admissible, whatever the variances
@test permutedims(t_cov, Val((2,1))) == transpose(t_cov)
@test permutedims(t_mix, Val((1,2))) === t_mix
@test permutedims(t_mix, Val((2,1))) == transpose(t_mix)
h3c = ThirdOrderTensorValue{cov2,cov2,3}(1:12...)
h3m = ThirdOrderTensorValue{con2,cov2,cov3}(1:12...)
@test permutedims(h3c, Val((1,2,3))) == h3c
@test permutedims(h3c, Val((2,1,3))) isa ThirdOrderTensorValue{cov2,cov2,3}
@test permutedims(h3c, Val((3,2,1))) isa ThirdOrderTensorValue{3,cov2,cov2}
@test permutedims(h3m, Val((3,1,2))) isa ThirdOrderTensorValue{cov3,con2,cov2}
let p = permutedims(h3m, Val((3,1,2)))
  @test all(p[c,a,b] == h3m[a,b,c] for a in 1:2, b in 1:2, c in 1:3)
end

# ==================================================================
# Symmetric and skew symmetric parts
# ==================================================================

@test symmetric_part(t_cov) isa SymTensorValue{cov2}
@test get_array(symmetric_part(t_cov)) == get_array(symmetric_part(t_int))
@test skew_symmetric_part(t_cov) isa SkewSymTensorValue{cov2}
@test_throws ErrorException symmetric_part(t_mix)
@test_throws ErrorException skew_symmetric_part(t_mix)
@test_throws ErrorException symmetric_part(TensorValue{cov2,2}(1,2,3,4))
# the specialised methods take a single index and are always admissible
@test symmetric_part(s_cov) === s_cov
@test skew_symmetric_part(k_cov) === k_cov
@test symmetric_part(k_cov) isa SymTensorValue{cov2,Int}
@test iszero(symmetric_part(k_cov))
@test skew_symmetric_part(s_cov) isa SkewSymTensorValue{cov2,Int}
@test iszero(skew_symmetric_part(s_cov))

# ==================================================================
# diag and broadcast
# ==================================================================

@test diag(t_int) == VectorValue(1,4)
@test_throws ErrorException diag(t_cov)
@test diag(k_int) == zero(VectorValue{2,Int})
@test_throws ErrorException diag(k_cov)

# Broadcasting pairs up components of the same index
@test (t_int .* t_int) == TensorValue{2,2}(1,4,9,16)
@test (t_cov .* t_cov) isa TensorValue{cov2,cov2}
@test (t_cov .* t_cov) == TensorValue{cov2,cov2}(1,4,9,16)
@test (s_cov .* s_cov) isa SymTensorValue{cov2}
@test_throws DimensionMismatch t_cov .* t_con
@test_throws DimensionMismatch t_cov .* t_int
@test_throws DimensionMismatch s_cov .* s_int

# Point with ContraVectorValue broadcasts
@test (v2_un .* v2) isa ContraVectorValue{con2}
@test (v2_un .* v2) == (v2 .* v2)
@test (v2 .* v2_un) == (v2 .* v2)
@test_throws DimensionMismatch v2_un .* ContraVectorValue{con3}((1.0,2.0,3.0))

# ==================================================================
# contracted_product
# ==================================================================

# The nine variance couples of a single contraction
for (Da, Db, ok) in ((2,2,true),      (2,cov2,false),   (2,con2,false),
                     (cov2,2,false),  (cov2,cov2,false),(cov2,con2,true),
                     (con2,2,false),  (con2,cov2,true), (con2,con2,false))
  a = TensorValue{cov2,Da}(1,2,3,4)
  b = TensorValue{Db,cov2}(1,2,3,4)
  if ok
    @test contracted_product(Val(1), a, b) isa TensorValue{cov2,cov2}
    @test get_array(contracted_product(Val(1), a, b)) ==
          get_array(contracted_product(Val(1), t_int, t_int))
  else
    @test_throws DimensionMismatch contracted_product(Val(1), a, b)
  end
end

# A zero-contraction never checks variance
@test (ω1 ⊗ ω1) isa TensorValue{cov2,cov2}
@test (ω1 ⊗ VectorValue(1,2)) isa TensorValue{cov2,2}
@test (VectorValue(1,2) ⊗ ω1) isa TensorValue{2,cov2}
@test (t_cov ⊗ t_cov) isa NTensorValue{Tuple{cov2,cov2,cov2,cov2}}

# Contracting away every index of one operand leaves the other's variance
@test (t_mix2 ⋅ ω1) isa CoVectorValue{cov2}
@test (ω13 ⋅ TensorValue{con3,cov3}(1:9...)) isa CoVectorValue{cov3}
# An order-1 contravariant result is an `ContraVectorValue`
@test (TensorValue{con2,con2}(1:4...) ⋅ ω1) isa ContraVectorValue{con2}
@test (TensorValue{con2,con2}(1:4...) ⋅ ω1) == ContraVectorValue{con2}((7,10))

# Dimension zero, contracted or not
z22 = TensorValue{cov2,con0,Int}(())
z02 = TensorValue{cov0,cov2,Int}(())
@test length(z22) == 0
@test contracted_product(Val(1), z22, z02) == zero(TensorValue{cov2,cov2,Int})
@test contracted_product(Val(1), z22, z02) isa TensorValue{cov2,cov2}
# the contracted zero dimensions must still be of opposite variance
@test_throws DimensionMismatch contracted_product(Val(1), TensorValue{cov2,cov0,Int}(()), z02)
# a zero dimension that is kept rather than contracted
zc2 = TensorValue{cov0,con2,Int}(())
@test contracted_product(Val(1), zc2, t_cov) isa TensorValue{cov0,cov2}
@test length(contracted_product(Val(1), zc2, t_cov)) == 0
# an empty operand still has its variance checked
@test_throws DimensionMismatch contracted_product(Val(1), z02, t_cov)

# ==================================================================
# tensor_contraction
# ==================================================================

# Two argument form, the nine couples on the contracted pair
for (Da, Db, ok) in ((2,2,true),      (2,cov2,false),   (2,con2,false),
                     (cov2,2,false),  (cov2,cov2,false),(cov2,con2,true),
                     (con2,2,false),  (con2,cov2,true), (con2,con2,false))
  a = TensorValue{Da,cov2}(1,2,3,4)
  b = TensorValue{cov2,Db}(1,2,3,4)
  if ok
    @test tensor_contraction(a, b, Val((1,)), Val((2,))) isa TensorValue{cov2,cov2}
  else
    @test_throws DimensionMismatch tensor_contraction(a, b, Val((1,)), Val((2,)))
  end
end

# Single argument form
@test tensor_contraction(t_mix, Val((1,)), Val((2,))) == tr(t_int)
@test tensor_contraction(t_mix2, Val((1,)), Val((2,))) == tr(t_int)
@test_throws DimensionMismatch tensor_contraction(t_cov, Val((1,)), Val((2,)))
@test_throws DimensionMismatch tensor_contraction(t_con, Val((1,)), Val((2,)))

h4 = NTensorValue{Tuple{con2,cov3,cov2,con3}}(reshape(1:36,2,3,2,3))
@test tensor_contraction(h4, Val((1,)), Val((3,))) isa TensorValue{cov3,con3}
@test tensor_contraction(h4, Val((2,)), Val((4,))) isa TensorValue{con2,cov2}
@test tensor_contraction(h4, Val((1,2)), Val((3,4))) isa Int
@test_throws DimensionMismatch tensor_contraction(h4, Val((1,)), Val((4,)))

# Kept indices keep their own variance, zero sized dimensions included
@test tensor_contraction(zc2, t_cov, Val((2,)), Val((1,))) isa TensorValue{cov0,cov2}
@test tensor_contraction(z22, z02, Val((2,)), Val((1,))) == zero(TensorValue{cov2,cov2,Int})

# ==================================================================
# SkewSymKTensorValue and its two variance aliases
# ==================================================================

# The aliases fix the variance parameter, so they are distinct types
@test ExteriorFormValue !== ExteriorVectorValue
@test CoVectorValue{3,Int,3} <: SkewSymKTensorValue
@test ContraVectorValue{3,Int,3} <: SkewSymKTensorValue
@test !(ContraVectorValue{3,Int,3} <: ExteriorFormValue)

# The constructors give the dimension the variance of the alias when it is
# given as a plain `Int`, and refuse to change an opposite one
ω_1 = CoVectorValue{3}((1,2,3))
v_1 = ContraVectorValue{3}((1,2,3))
@test typeof(ω_1) === CoVectorValue{cov3,Int,3}
@test typeof(v_1) === ContraVectorValue{con3,Int,3}
@test typeof(ContraVectorValue{con3}((1,2,3))) === typeof(v_1)
@test_throws AssertionError ContraVectorValue{cov3}((1,2,3))
@test typeof(SkewSymKTensorValue{Con,1,3,Int,3}(1,2,3)) === typeof(v_1)
@test size(v_1) == (con3,)
@test tensor_type(typeof(v_1)) == (1,0,0)
@test tensor_type(typeof(ExteriorVectorValue{2,3}((1,2,3)))) == (2,0,0)

@test_throws AssertionError ContraVectorValue{3}((1,2))
@test_throws AssertionError SkewSymKTensorValue{Con,1,3,Int,3}(1,2)

# The whole interface follows the variance
@test num_indep_components(ExteriorVectorValue{2,con3}) == 3
@test zero(ExteriorVectorValue{2,3,Float64}) === ExteriorVectorValue{2,con3}((0.0,0.0,0.0))
@test one(ExteriorVectorValue{0,3,Float64}) === 1.0
@test change_eltype(typeof(v_1), Float64) === ContraVectorValue{con3,Float64,3}
@test change_eltype(typeof(ω_1), Float64) === CoVectorValue{cov3,Float64,3}

# Indexing and expansion to the full tensor do not depend on the variance
v_2 = ExteriorVectorValue{2,3}((1,2,3))
@test (v_2[1,2], v_2[2,1], v_2[1,1]) == (1, -1, 0)
@test get_array(v_2) == [0 1 2; -1 0 3; -2 -3 0]
@test convert(SArray{Tuple{3,3},Int}, v_2) == get_array(v_2)

# A form and a vector of the same order contract through `⋅`; `inner` wants the
# two in one tensor space and so refuses them, varianced pair or not
@test (v_1 ⋅ ω_1) == 14
@test outer(v_1, ω_1) isa TensorValue{con3,cov3}
@test_throws ErrorException inner(v_1, ω_1)
@test_throws ErrorException inner(v_2, ExteriorFormValue{2,3}((1,2,3)))
@test_throws ArgumentError inner(v_1, v_1)

# Values of the two variances mix in no other way
@test_throws ErrorException v_1 + ω_1
@test_throws ErrorException convert(CoVectorValue{3,Float64}, v_1)

# A vector is displayed on the frame dual to the coframe of a form
_str_v(io_props, a) = sprint((io, x) -> show(io, MIME("text/plain"), x), a;
                             context = io_props)
@test _str_v(:coordinates => :cartesian,   v_1) == "1 ∂₁ + 2 ∂₂ + 3 ∂₃"
@test _str_v(:coordinates => :barycentric, v_1) == "1 ∂λ₁ + 2 ∂λ₂ + 3 ∂λ₃"
@test _str_v(:coordinates => :cartesian,   ω_1) == "1 dx¹ + 2 dx² + 3 dx³"
@test occursin("∂₁ ∧ ∂₂", _str_v(:coordinates => :cartesian, v_2))
@test_throws ArgumentError _str_v(:coordinates => :something_else, v_1)

@test indep_components_names(ContraVectorValue{con3}) == ["∂₁", "∂₂", "∂₃"]
@test indep_components_names(ExteriorVectorValue{2,con3}) ==
      ["∂₁ ∧ ∂₂", "∂₁ ∧ ∂₃", "∂₂ ∧ ∂₃"]

end # module
