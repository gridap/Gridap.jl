module HigherOrderGradientsTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Fields: ad_return_cache, ad_evaluate!
using Gridap.Polynomials

# PolynomialBasis: hand-written kernels for N ≤ 2, forward-mode AD above. The AD evaluation
# reproduces the hand-written kernels exactly where both exist.
for (D, k) in ((2, 3), (3, 2))
  b = MonomialBasis(Val(D), Float64, k)
  xs = [Point(ntuple(d -> 0.1d + 0.05i, D)) for i in 1:4]
  for N in 1:2
    fg = FieldGradientArray{N}(b)
    @test ad_evaluate!(ad_return_cache(fg, xs), fg, xs) == evaluate(fg, xs)
  end
  for N in 3:4
    fg = FieldGradientArray{N}(b)
    r = evaluate(fg, xs)
    @test size(r) == (length(xs), length(b))
    @test eltype(r) == Fields.nth_gradient_type(Float64, first(xs), Val(N))
    # The N-th gradient is the gradient of the (N-1)-th
    rp = evaluate(Broadcasting(∇)(FieldGradientArray{N-1}(b)), xs)
    @test r ≈ rp
    # Single point
    @test evaluate(fg, first(xs)) ≈ r[1,:]
  end
end

end # module
