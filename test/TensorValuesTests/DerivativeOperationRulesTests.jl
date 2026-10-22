module DerivativeOperationRulesTests

using Gridap.TensorValues
using Gridap.TensorValues: _PartialsTensor, _map_second_factor, _leibniz_cross_term
using Gridap.Arrays: get_array
using Test

# mock partials
∂(d, b) = product_tensor(d, b)
as_factor(b::MultiValue) = b
as_factor(b::Number) = ScalarValue(b)

# Component array of a value, whatever its tensor order, including none
full(x::MultiValue) = get_array(x)
full(x::Number) = fill(x)

const d1 = SymNTensorValue{Tuple{Cov(2)},Float64}(1.0, 0.5)
const e1 = SymNTensorValue{Tuple{Cov(2)},Float64}(0.25, 2.0)
const d2 = SymNTensorValue{Tuple{Cov(2),Cov(2)},Float64}(1.0, 0.5, 2.0)
const e2 = SymNTensorValue{Tuple{Cov(2),Cov(2)},Float64}(0.5, 1.5, 0.25)

# Bilinear operations with operand values and jet second factors they accept
const BILINEAR_CASES = (
  (*, ScalarValue(2.0), VectorValue(3.0,4.0),
      ScalarValue(0.3), VectorValue(0.1,0.2)),
  (⋅, VectorValue(1.0,2.0), VectorValue(3.0,4.0),
      VectorValue(0.1,0.2), VectorValue(0.3,0.1)),
  (⊙, TensorValue(1.0,2.0,3.0,4.0), TensorValue(2.0,1.0,0.5,1.0),
      TensorValue(0.1,0.2,0.1,0.3), TensorValue(0.2,0.1,0.3,0.1)),
  (⊗, VectorValue(1.0,2.0), VectorValue(3.0,4.0),
      VectorValue(0.1,0.2), VectorValue(0.3,0.1)),
  (∧, ExteriorFormValue{1,2,Float64}(1.0,2.0), ExteriorFormValue{1,2,Float64}(3.0,4.0),
      ExteriorFormValue{1,2,Float64}(0.1,0.2), ExteriorFormValue{1,2,Float64}(0.3,0.1)),
  (ι, ContraVectorValue{Con(2)}(1.0,2.0), ExteriorFormValue{2,Cov(2),Float64}(3.0),
      ContraVectorValue{Con(2)}(0.1,0.2), ExteriorFormValue{2,Cov(2),Float64}(0.5)),
  (koszul, ContraVectorValue{Con(2)}(1.0,2.0), ExteriorFormValue{2,Cov(2),Float64}(3.0),
      ContraVectorValue{Con(2)}(0.1,0.2), ExteriorFormValue{2,Cov(2),Float64}(0.5)),
  (alt_endomorphism, TensorValue(2.0,1.0,0.5,3.0), ExteriorFormValue{1,Cov(2),Float64}(1.0,2.0),
      TensorValue(0.1,0.2,0.1,0.3), ExteriorFormValue{1,Cov(2),Float64}(0.3,0.1)),
)

# The first order expansion
#   op(u+δu, v+δv) - op(u,v) = op(δu,v) + op(u,δv) + op(δu,δv)
# validates the order 1 rule: taking `δ` of magnitude `h` along each derivation
# direction, the neglected term is of magnitude `h²`, so the finite difference
# matches the rule to `h`.
function first_order_discrepancy(op, u, v, up, vp, du, dv; h=1.0e-3)
  R = get_array(partials_op_rule(Val(1),op)(u, v, ∂(du,up), ∂(dv,vp)))
  worst = 0.0
  for i in axes(R,1)
    δu, δv = (h*du[i])*up, (h*dv[i])*vp
    fd = (full(op(u+δu, v+δv)) .- full(op(u,v))) ./ h
    worst = max(worst, maximum(abs.(selectdim(R,1,i) .- fd)))
  end
  worst
end

# The residual the first order expansion leaves is the mixed term itself,
# `op(δu,δv)`, which the order 2 rule holds at `(i,j)` as
# `op(∂ᵢu,∂ⱼv) + op(∂ⱼu,∂ᵢv)`. The operation being bilinear, the only error is
# the cancellation of the `h²` scaling.
function cross_term_discrepancy(op, u, v, up, vp, du, dv; h=1.0e-2)
  C = get_array(_leibniz_cross_term(op, ∂(du,up), ∂(dv,vp)))
  worst = 0.0
  for i in axes(C,1), j in axes(C,2)
    residual(a, b) = full(op(u+a, v+b)) .- full(op(u,v)) .- full(op(a,v)) .- full(op(u,b))
    r = residual((h*du[i])*up, (h*dv[j])*vp) .+ residual((h*du[j])*up, (h*dv[i])*vp)
    worst = max(worst, maximum(abs.(selectdim(selectdim(C,1,i),1,j) .- r ./ h^2)))
  end
  worst
end

# The alias qualifying partial derivatives

# It reads back the derivation order, the coordinate dimension and the value type
@test ∂(d1, VectorValue(1.0,2.0)) isa _PartialsTensor{1,Cov(2),VectorValue{2}}
@test ∂(d2, VectorValue(1.0,2.0)) isa _PartialsTensor{2,Cov(2),VectorValue{2}}
@test !(∂(d2, VectorValue(1.0,2.0)) isa _PartialsTensor{1,Cov(2)})

# A product tensor whose first factor is not a derivation block holds no partials
@test !(product_tensor(VectorValue(1.0,2.0), VectorValue(1.0,2.0)) isa _PartialsTensor)

# Jointly linear operations: the rule is the operation itself, at any order

@test partials_op_rule(+) === +
@test partials_op_rule(-) === -

# The sum of two jets along the same directions differentiates their sum
@test partials_op_rule(+)(∂(d1,VectorValue(1.0,2.0)), ∂(d1,VectorValue(0.5,1.0))) ===
      ∂(d1, VectorValue(1.5,3.0))
@test partials_op_rule(-)(∂(d2,VectorValue(1.0,2.0)), ∂(d2,VectorValue(0.5,1.0))) ===
      ∂(d2, VectorValue(0.5,1.0))

# ⋆ is linear with constant coefficients, so it acts on the value factor alone,
# at any derivation order
let ω = ExteriorFormValue{1,2,Float64}(1.0,2.0), ⋆ω = flat_hodge_star(ω)
  @test partials_op_rule(flat_hodge_star)(∂(d1,ω)) === ∂(d1, ⋆ω)
  @test partials_op_rule(flat_hodge_star)(∂(d2,ω)) === ∂(d2, ⋆ω)
  @test partials_op_rule(flat_hodge_star)(∂(d1,ω)) === _map_second_factor(flat_hodge_star, ∂(d1,ω))
end

# Bilinear operations, order 1

for (op, u, v, up, vp) in BILINEAR_CASES
  rule = partials_op_rule(Val(1), op)

  # Each argument's jet acts against the other argument's value, in its own
  # position: a swapped argument changes the result of a non commutative `op`
  expected = ∂(d1, as_factor(op(up,v))) + ∂(e1, as_factor(op(u,vp)))
  @test typeof(rule(u, v, ∂(d1,up), ∂(e1,vp))) === typeof(expected)
  @test rule(u, v, ∂(d1,up), ∂(e1,vp)) ≈ expected

  @test first_order_discrepancy(op, u, v, up, vp, d1, e1) < 1.0e-2
end

# A scalar field's value is a plain number, where its own derivatives hold a
# `ScalarValue`, and the rule accepts either
let v = VectorValue(3.0,4.0), ∂u = ∂(d1,ScalarValue(0.3)), ∂v = ∂(e1,VectorValue(0.1,0.2))
  @test partials_op_rule(Val(1),*)(2.0, v, ∂u, ∂v) ===
        partials_op_rule(Val(1),*)(ScalarValue(2.0), v, ∂u, ∂v)
end

# Bilinear operations, order 2

for (op, u, v, up, vp) in BILINEAR_CASES
  rule = partials_op_rule(Val(2), op)
  ∂u, ∂v = ∂(d1,up), ∂(e1,vp)
  ∂∂u, ∂∂v = ∂(d2,up), ∂(e2,vp)

  # The three terms are independent: vanishing second derivatives leave the
  # mixed term alone
  @test rule(u, v, ∂u, ∂v, zero(∂∂u), zero(∂∂v)) === _leibniz_cross_term(op, ∂u, ∂v)

  # and vanishing first derivatives leave the two outer terms, each second
  # derivative acting against the other argument's value
  outer_terms = ∂(d2, as_factor(op(up,v))) + ∂(e2, as_factor(op(u,vp)))
  @test typeof(rule(u, v, zero(∂u), zero(∂v), ∂∂u, ∂∂v)) === typeof(outer_terms)
  @test rule(u, v, zero(∂u), zero(∂v), ∂∂u, ∂∂v) ≈ outer_terms

  @test cross_term_discrepancy(op, u, v, up, vp, d1, e1) < 1.0e-2
end

# The second factor action

# The second factor at a given independent component of the first is that component
# times the second factor, so a linear action carries through the whole product
av, bs = VectorValue(1.0,2.0), SymTensorValue(1.0,2.0,3.0)
pm = product_tensor(av,bs)
@test _map_second_factor(identity, pm) === pm
@test _map_second_factor(b -> 2*b, pm) === 2*pm

# A scalar valued action gives a ScalarValue second factor, so the result is a
# product tensor whatever the action returns
@test _map_second_factor(tr, pm) === product_tensor(av, ScalarValue(tr(bs)))

# The second factor type follows the action
pt = product_tensor(av, TensorValue(1.0,2.0,3.0,4.0))
@test _map_second_factor(symmetric_part, pt) === product_tensor(av, symmetric_part(TensorValue(1.0,2.0,3.0,4.0)))

# Independent components are read layout agnostically, so a factor storing dependent
# components is mapped like any other, in either position
bst = SymTracelessTensorValue(1.0,2.0)
@test _map_second_factor(b -> 2*b, product_tensor(av,bst)) === product_tensor(av,2*bst)
@test _map_second_factor(b -> 2*b, product_tensor(bst,av)) === product_tensor(bst,2*av)
@test _map_second_factor(identity, product_tensor(bst,bst)) === product_tensor(bst,bst)

# The action reaches an order 0 second factor like any other
p0 = product_tensor(VectorValue(1.0,2.0), ScalarValue(3.0))
@test _map_second_factor(b -> 2*b, p0) === 2*p0

# A factor with no component gives the action no argument
@test_throws ErrorException _map_second_factor(identity, product_tensor(zero(VectorValue{0,Float64}), bs))

end # module
