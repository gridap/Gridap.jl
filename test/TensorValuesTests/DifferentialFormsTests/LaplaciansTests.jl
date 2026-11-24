module LaplaciansTests

using Test
using Gridap.TensorValues
import LinearAlgebra

# Deterministic components, of either sign
w(n) = sin(1.7n) + 0.3n - 0.5

# Mock second order partial derivatives of `a`, shaped as
# `second_partial_derivatives` produces: the two covariant derivation indices
# share the symmetric first factor. The identities below follow from that
# symmetry and from the linearity of `𝑑` and `δ` alone, so an arbitrary `∂∂a`
# exercises them.
function ∂∂mock(V, D, seed=0)
  A = SymNTensorValue{Tuple{D,D}}
  L = num_indep_components(V{Float64}) * num_indep_components(A{Float64})
  ProductTensorValue{A,V{Float64},Float64}(ntuple(n -> w(n+seed), L))
end

# `(dδω)_I = -∑_q (-1)^(q-1) ∑_j ∂∂ω[I_q, j, j, I∖I_q]`, `δω` contracting the
# derivation index of `∂ω` with the first index of `ω`, and `𝑑` alternating.
function ref_dδ(∂∂ω, D, K, I)
  sum(1:K; init=0.0) do q
    Iq = ntuple(i -> I[i + Int(i≥q)], K-1)
    -(-1)^(q-1) * sum(j -> ∂∂ω[I[q], j, j, Iq...], 1:D)
  end
end

# `(δdω)_I = -∑_j ∑_a (-1)^(a-1) ∂∂ω[j, M_a, M∖M_a]` with `M = (j, I...)`, the
# degree raised by `𝑑` before the contraction.
function ref_δd(∂∂ω, D, K, I)
  sum(1:D; init=0.0) do j
    M = (j, I...)
    sum(1:K+1; init=0.0) do a
      Ma = ntuple(i -> M[i + Int(i≥a)], K)
      -(-1)^(a-1) * ∂∂ω[j, M[a], Ma...]
    end
  end
end

ref_bochner(∂∂a, D, I) = -sum(j -> ∂∂a[j, j, I...], 1:D)

# `ΔdR` is `dδ + δd` and `ΔB` is the componentwise laplacian, degree by degree
# and index by index of the full antisymmetric tensors.
function test_form_laplacians(K, D)
  Dc = Cov(D)
  ∂∂ω = ∂∂mock(ExteriorFormValue{K,Dc}, Dc, K+D)
  ΔdRω = second_partials_to_flat_hodgelaplacian(∂∂ω)
  ΔBω  = second_partials_to_flat_bochnerlaplacian(∂∂ω)

  @test ΔdRω === ΔBω
  @test ΔdRω isa ExteriorFormValue{K,Dc,Float64}
  @test ΔdRω === @inferred second_partials_to_flat_hodgelaplacian(∂∂ω)

  for I in Iterators.product(ntuple(_ -> 1:D, K)...)
    @test ΔdRω[I...] ≈ ref_dδ(∂∂ω, D, K, I) + ref_δd(∂∂ω, D, K, I) atol=1.0e-13
    @test ΔdRω[I...] ≈ ref_bochner(∂∂ω, D, I) atol=1.0e-13
  end
end

for D in (2,3), K in 0:D
  test_form_laplacians(K, D)
end

# A scalar field is the 0-form whose `δ` has no value type, so `ΔdR` is `δ𝑑`
let ∂∂f = ∂∂mock(ScalarValue, Cov(3), 5)
  @test second_partials_to_flat_hodgelaplacian(∂∂f) === second_partials_to_flat_bochnerlaplacian(∂∂f)
  @test second_partials_to_flat_bochnerlaplacian(∂∂f) ≈ ref_bochner(∂∂f, 3, ())
  @test second_partials_to_flat_bochnerlaplacian(∂∂f) isa Float64
end

# `ΔB` takes a tensor of any order and variance, and keeps its symmetries
let ∂∂s = ∂∂mock(SymTensorValue{Cov(3)}, Cov(3), 2),
    ∂∂t = ∂∂mock(TensorValue{Cov(2),Con(2)}, Cov(2), 4)
  @test second_partials_to_flat_bochnerlaplacian(∂∂s) isa SymTensorValue{Cov(3),Float64}
  @test second_partials_to_flat_bochnerlaplacian(∂∂t) isa TensorValue{Cov(2),Con(2),Float64}
  for i in 1:3, j in 1:3
    @test second_partials_to_flat_bochnerlaplacian(∂∂s)[i,j] ≈ ref_bochner(∂∂s, 3, (i,j))
  end
  for i in 1:2, j in 1:2
    @test second_partials_to_flat_bochnerlaplacian(∂∂t)[i,j] ≈ ref_bochner(∂∂t, 2, (i,j))
  end
end

# The metric aware kernels take their metric as `δ` does. Fed the euclidean one,
# whose Christoffel symbols vanish with the partial derivatives of `g`, they
# answer what the flat kernels above answer. Their curved values are checked
# against `dδ + δd` of a varying metric in `FieldsTests`.
function ∂mock(V, D, seed=0)
  A = SymNTensorValue{Tuple{D}}
  L = num_indep_components(V{Float64}) * num_indep_components(A{Float64})
  ProductTensorValue{A,V{Float64},Float64}(ntuple(n -> w(n+seed), L))
end

function test_flat_metric(K, D)
  Dc = Cov(D)
  V = ExteriorFormValue{K,Dc}
  ∂∂ω = ∂∂mock(V, Dc, K+D)
  ∂ω  = ∂mock(V, Dc, K+D+5)
  ω   = V{Float64}(ntuple(n -> w(n+2), binomial(D,K)))
  ∂g  = zero(ProductTensorValue{SymNTensorValue{Tuple{Dc}},SymTensorValue{Dc},Float64})
  ∂∂g = zero(ProductTensorValue{SymNTensorValue{Tuple{Dc,Dc}},SymTensorValue{Dc},Float64})
  args = (∂∂ω, ∂ω, ω, ∂∂g, ∂g, LinearAlgebra.I)

  ΔBω = second_partials_to_bochnerlaplacian(args...)
  @test ΔBω isa ExteriorFormValue{K,Dc,Float64}
  @test ΔBω ≈ second_partials_to_flat_bochnerlaplacian(∂∂ω)
  @test second_partials_to_hodgelaplacian(args...) ≈ second_partials_to_flat_hodgelaplacian(∂∂ω)
end

for (K,D) in ((1,2), (2,3))
  test_flat_metric(K,D)
end

# Partial derivatives are covariant, a contravariant derivation index is not one
@test_throws ArgumentError second_partials_to_flat_bochnerlaplacian(
  ∂∂mock(ExteriorFormValue{1,Cov(2)}, Con(2)))

# `ΔdR` is the laplacian of a form, a tensor of two indices is not one
@test_throws MethodError second_partials_to_flat_hodgelaplacian(
  ∂∂mock(SymTensorValue{Cov(2)}, Cov(2)))

end # module
