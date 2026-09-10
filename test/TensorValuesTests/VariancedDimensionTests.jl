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

ω13 = ExteriorFormValue{1,3}((1,2,3))
@test size(ω13) == (cov3,)
@test ω13 isa ExteriorFormValue{1,cov3}
@test ω13 ⊗ ω13 == TensorValue{cov3, cov3}(1, 2, 3, 2, 4, 6, 3, 6, 9)

ω23 = ExteriorFormValue{2,3}((1,2,3))
@test size(ω23) == (cov3,cov3)

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

ω1      = ExteriorFormValue{1,cov2}((1,2))

# ==================================================================
# Comparison
# ==================================================================

@test t_mix ≈ t_mix
@test t_int ≈ TensorValue{2,2}(1.0,2.0,3.0,4.0)
# `get_array` erases variance: sizes equal in value but not in variance compare false
@test !(t_mix ≈ t_int)
@test !(t_mix ≈ t_cov)
@test !(t_mix ≈ t_mix2)
@test !(s_cov ≈ s_int)
# a genuine size mismatch is still reported by the array comparison
@test_throws DimensionMismatch (t_int ≈ TensorValue{2,3}(1:6...))

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
@test_throws DimensionMismatch g \ ExteriorFormValue{1,cov2}((1,2))
# left operand not square
@test_throws DimensionMismatch TensorValue{2,3}(1:6...) \ VectorValue(1,2)
# an order-0 right operand is not supported
@test_throws ErrorException t_int \ ExteriorFormValue{0,2}((1,))

# `inv` of a (0,2) tensor is a (2,0) tensor, so dividing a covariant vector by it
# would give a contravariant vector, which Gridap has no type for
@test_throws ErrorException t_cov \ ExteriorFormValue{1,cov2}((1,2))
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

# Symmetric tensors
@test s_int ⊙ s_int == 1 + 2*(2*2) + 3*3
@test s_cov ⊙ s_con == s_int ⊙ s_int
@test s_con ⊙ s_cov == s_int ⊙ s_int
@test_throws DimensionMismatch s_cov ⊙ s_cov
@test_throws DimensionMismatch s_con ⊙ s_con
@test_throws DimensionMismatch s_cov ⊙ SymTensorValue{con3}(1:6...)

# Traceless tensors go through the same methods
@test q_cov ⊙ q_con == q_int ⊙ q_int
@test_throws DimensionMismatch q_cov ⊙ q_cov

# Skew symmetric tensors
@test k_cov ⊙ k_con == k_int ⊙ k_int
@test_throws DimensionMismatch k_cov ⊙ k_cov

# A skew symmetric tensor is orthogonal to every symmetric one
@test iszero(k_int ⊙ s_int)
@test iszero(k_cov ⊙ s_con)
@test iszero(s_con ⊙ k_cov)
@test_throws DimensionMismatch k_cov ⊙ s_cov
@test_throws DimensionMismatch s_cov ⊙ k_cov

# Fourth order symmetric tensors
@test f_cov ⊙ f_con == f_int ⊙ f_int
@test_throws DimensionMismatch f_cov ⊙ f_cov

# The generic method pairs the indices of equal rank, whatever the sizes
@test_throws DimensionMismatch t_cov ⊙ t_cov
@test t_mix ⊙ t_mix2 == t_int ⊙ t_int
@test contracted_product(Val(2), t_mix, t_mix2) == t_int ⊙ t_int
# operands of different order are still refused before any index is paired
@test_throws ErrorException ω1 ⊙ t_cov
# same order but incompatible dimensions is a contraction failure
@test_throws DimensionMismatch VectorValue(1,2,3) ⊙ VectorValue(1,2)

# The order-4 against order-2 methods defer their variance check downstream
@test get_array(f_cov ⊙ s_con) == get_array(f_int ⊙ s_int)
@test get_array(s_con ⊙ f_cov) == get_array(s_int ⊙ f_int)
@test_throws DimensionMismatch f_cov ⊙ s_cov
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
@test (f_cov ⋅² f_con) isa HighOrderTensorValue{Tuple{cov2,cov2,con2,con2}}
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
@test tr(h3) isa ExteriorFormValue{1,cov3}
@test get_array(tr(h3)) == get_array(tr(h3i))
@test_throws ArgumentError tr(ThirdOrderTensorValue{cov2,cov2,3}(1:12...))

# `tr` and the self-contraction agree, variance included
@test tr(h3) == tensor_contraction(h3, Val(1), Val(2))
@test tr(t_mix) == tensor_contraction(t_mix, Val(1), Val(2))
@test_throws DimensionMismatch tensor_contraction(t_cov, Val(1), Val(2))

# ==================================================================
# meas and norm
# ==================================================================

@test meas(t_mix) == meas(t_int)
@test_throws ErrorException meas(t_cov)
# a row read as a vector is only length preserving without variance
@test meas(TensorValue{1,2}(3,4)) == 5
@test_throws ErrorException meas(TensorValue{1,cov2}(3,4))

# `norm` inherits its variance rules from `inner`, which pairs index k of `u`
# with index k of `u`: no varianced tensor has an invariant Frobenius norm, a
# (1,1) tensor included -- the invariant quantity there is `tr(u⋅u)`.
@test norm(t_int) == sqrt(1 + 4 + 9 + 16)
@test_throws DimensionMismatch norm(t_cov)
@test_throws DimensionMismatch norm(t_mix)
@test tr(t_mix ⋅ t_mix) == tr(t_int ⋅ t_int)
# a form carries covariant indices, so its Euclidean norm goes through
# `flat_inner`, which raises them by the identity metric
@test_throws DimensionMismatch norm(ω1)
@test sqrt(flat_inner(ω1, ω1)) == norm(VectorValue(1,2))
@test ω1 ⊙ᴱ ω1 === flat_inner(ω1, ω1)

# ==================================================================
# Variance of a metric tensor
# ==================================================================

# A metric lowers indices and is a (0,2) tensor; its inverse raises them and is
# a (2,0) one. An un-varianced metric is accepted, as everywhere else.
ω2   = ExteriorFormValue{1,cov2}((1.0, 2.0))
v2   = VectorValue(1.0, 2.0)
g_id = SymTensorValue{2}(1.0, 0.0, 1.0)
g_lo = SymTensorValue{cov2}(1.0, 0.0, 1.0)
g_up = SymTensorValue{con2}(1.0, 0.0, 1.0)

@test tensor_type(g_lo) == (0,2,0)
@test tensor_type(g_up) == (2,0,0)
@test tensor_type(g_id) == (0,0,2)

# `hodge_star`, `sharp` and `form_inner` raise indices: they take `g_inv`
@test hodge_star(ω2, g_up, 1.0) ≈ flat_hodge_star(ω2)
@test sharp(ω2, g_up) ≈ v2
@test form_inner(ω2, ω2, g_up) ≈ flat_form_inner(ω2, ω2)
for g_ko in (g_lo, g_id)
  @test_throws ErrorException hodge_star(ω2, g_ko, 1.0)
  @test_throws ErrorException sharp(ω2, g_ko)
  @test_throws ErrorException form_inner(ω2, ω2, g_ko)
end

# `flat` lowers an index: it takes `g`
@test flat(v2, g_lo) ≈ ω2
@test_throws ErrorException flat(v2, g_id)
@test_throws ErrorException flat(v2, g_up)

# metric errors
@test_throws "(2, 0)" sharp(ω2, g_lo)
@test_throws "(0, 2)" flat(v2, g_up)
@test_throws DimensionMismatch sharp(ω2, SymTensorValue{con3}(1.,0.,0.,1.,0.,1.))
@test_throws DimensionMismatch flat(v2, SymTensorValue{cov3}(1.,0.,0.,1.,0.,1.))

# The result of `flat` carries the covariant dimension of the metric
@test flat(v2, g_lo) isa ExteriorFormValue{1,cov2}
@test hodge_star(ω2, g_up, 1.0) isa ExteriorFormValue{1,cov2}
@test sharp(ω2, g_up) isa VectorValue{2}

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
@test (t_cov ⊗ t_cov) isa HighOrderTensorValue{Tuple{cov2,cov2,cov2,cov2}}

# Contracting away every index of one operand leaves the other's variance
@test (t_mix2 ⋅ ω1) isa ExteriorFormValue{1,cov2}
@test (ω13 ⋅ TensorValue{con3,cov3}(1:9...)) isa ExteriorFormValue{1,cov3}
# Gridap has no contravariant vector, so an order-1 contravariant result throws
@test_throws ErrorException TensorValue{con2,con2}(1:4...) ⋅ ω1

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

h4 = HighOrderTensorValue{Tuple{con2,cov3,cov2,con3}}(reshape(1:36,2,3,2,3))
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
@test ExteriorFormValue{1,3,Int,3} <: SkewSymKTensorValue
@test ExteriorVectorValue{1,3,Int,3} <: SkewSymKTensorValue
@test !(ExteriorVectorValue{1,3,Int,3} <: ExteriorFormValue)

# The constructors give the dimension the variance of the alias when it is
# given as a plain `Int`, and refuse to change an opposite one
ω_1 = ExteriorFormValue{1,3}((1,2,3))
v_1 = ExteriorVectorValue{1,3}((1,2,3))
@test typeof(ω_1) === ExteriorFormValue{1,cov3,Int,3}
@test typeof(v_1) === ExteriorVectorValue{1,con3,Int,3}
@test typeof(ExteriorVectorValue{1,con3}((1,2,3))) === typeof(v_1)
@test_throws AssertionError ExteriorVectorValue{1,cov3}((1,2,3))
@test typeof(SkewSymKTensorValue{1,3,Int,3,Con}(1,2,3)) === typeof(v_1)
@test size(v_1) == (con3,)
@test tensor_type(typeof(v_1)) == (1,0,0)
@test tensor_type(typeof(ExteriorVectorValue{2,3}((1,2,3)))) == (2,0,0)

@test_throws AssertionError ExteriorVectorValue{1,3}((1,2))
@test_throws AssertionError SkewSymKTensorValue{1,3,Int,3,Con}(1,2)

# The whole interface follows the variance
@test num_indep_components(ExteriorVectorValue{2,con3}) == 3
@test zero(ExteriorVectorValue{2,3,Float64}) === ExteriorVectorValue{2,con3}((0.0,0.0,0.0))
@test one(ExteriorVectorValue{0,3,Float64}) === 1.0
@test change_eltype(typeof(v_1), Float64) === ExteriorVectorValue{1,con3,Float64,3}
@test change_eltype(typeof(ω_1), Float64) === ExteriorFormValue{1,cov3,Float64,3}

# Indexing and expansion to the full tensor do not depend on the variance
v_2 = ExteriorVectorValue{2,3}((1,2,3))
@test (v_2[1,2], v_2[2,1], v_2[1,1]) == (1, -1, 0)
@test get_array(v_2) == [0 1 2; -1 0 3; -2 -3 0]
@test convert(SArray{Tuple{3,3},Int}, v_2) == get_array(v_2)

# A form and a vector of the same order contract, two forms do not
@test inner(v_1, ω_1) == 14
@test inner(v_2, ExteriorFormValue{2,3}((1,2,3))) == 28
@test (v_1 ⋅ ω_1) == 14
@test outer(v_1, ω_1) isa TensorValue{con3,cov3}
@test_throws DimensionMismatch inner(v_1, v_1)

# Values of the two variances mix in no other way
@test_throws ErrorException v_1 + ω_1
@test_throws ErrorException convert(ExteriorFormValue{1,3,Float64}, v_1)
# the pointwise algebra is defined on forms only
@test_throws MethodError v_1 ∧ v_1
@test_throws MethodError flat_hodge_star(v_1)

# A vector is displayed on the frame dual to the coframe of a form
_str_v(io_props, a) = sprint((io, x) -> show(io, MIME("text/plain"), x), a;
                             context = io_props)
@test _str_v(:coordinates => :cartesian,   v_1) == "1 ∂₁ + 2 ∂₂ + 3 ∂₃"
@test _str_v(:coordinates => :barycentric, v_1) == "1 ∂λ₁ + 2 ∂λ₂ + 3 ∂λ₃"
@test _str_v(:coordinates => :cartesian,   ω_1) == "1 dx¹ + 2 dx² + 3 dx³"
@test occursin("∂₁ ∧ ∂₂", _str_v(:coordinates => :cartesian, v_2))
@test_throws ArgumentError _str_v(:coordinates => :something_else, v_1)

@test indep_components_names(ExteriorVectorValue{1,con3}) == ["∂₁", "∂₂", "∂₃"]
@test indep_components_names(ExteriorVectorValue{2,con3}) ==
      ["∂₁ ∧ ∂₂", "∂₁ ∧ ∂₃", "∂₂ ∧ ∂₃"]

end # module
