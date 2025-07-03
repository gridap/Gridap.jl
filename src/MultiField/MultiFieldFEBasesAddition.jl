# This code is written solely to use `MultiField` and keep the form of a sum of spaces.
struct MultiFieldFEBasesAddition{A, B}
	a::A
	b::B
end

function Base.:+(a::MultiFieldFEBasisComponent, b::MultiFieldFEBasisComponent)
	return MultiFieldFEBasesAddition(a, b)
end

for op in (:inner, :dot, :*, :outer)
	@eval TensorValues.$op(mfa::MultiFieldFEBasesAddition, f) = $op(f, mfa)
	@eval TensorValues.$op(f, mfa::MultiFieldFEBasesAddition) = $op(f, mfa.a) + $op(f, mfa.b)
	@eval TensorValues.$op(mfa::MultiFieldFEBasesAddition, b::MultiFieldFEBasesAddition) = $op(mfa.a, b) + $op(mfa.b, b)
end
for op in (:gradient, :symmetric_gradient)
	@eval Fields.$op(mfa::MultiFieldFEBasesAddition) = $op(mfa.a) + $op(mfa.b)
end
for op in (:adjoint,)
	@eval Base.$op(mfa::MultiFieldFEBasesAddition) = $op(mfa.a) + $op(mfa.b)
end
