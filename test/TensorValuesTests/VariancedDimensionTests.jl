module VariancedDimensionTests

using Test
using Gridap
using Gridap.TensorValues

t = TensorValue{2, cov2}(1, 2, 3, 4)
@test 2 * t == TensorValue{2, Covariant(2)}(2, 4, 6, 8)
@test t + t == 2 * t

t2 = TensorValue{2, 2}(1:4...)
@test_throws "failed to change any argument" t + t2

v = VectorValue(5,7)
v_outer_t = ThirdOrderTensorValue{2, 2, cov2}(5, 7, 10, 14, 15, 21, 20, 28)
@test v ⊗ t == v_outer_t
@test t ⊗ v != v_outer_t

ω13 = ExteriorFormValue{1,3}((1,2,3))
@test size(ω13) == (cov3,)
@test ω13 isa ExteriorFormValue{1,cov3}
@test ω13 ⊗ ω13 == TensorValue{cov3, cov3}(1, 2, 3, 2, 4, 6, 3, 6, 9)

ω23 = ExteriorFormValue{2,3}((1,2,3))
@test size(ω23) == (cov3,cov3)

@test tensor_type(1) == (0,0)
@test tensor_type(v) == (1,0)
@test tensor_type(t) == (1,1)
@test tensor_type(ω13) == (0,1)
@test tensor_type(ω23) == (0,2)
@test tensor_type(v_outer_t) == (2,1)

end # module
