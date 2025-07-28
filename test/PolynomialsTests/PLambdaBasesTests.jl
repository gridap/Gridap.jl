module PLambdaBasisTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using ForwardDiff
using StaticArrays

using Gridap.Polynomials: _combination_index, bernstein_term_id

r = 4

# Bubble indices validation

D = 3
N = D+1
for k in 0:D
  w_prev = 0
  for (F, bubble_functions) in PmΛ_bubbles(r,k,D)
    @test issorted(F)  &&  (k ≤ length(F) ≤ N)  &&  (F ⊆ 1:N) || (k,D,F)

    passed = true
    for (w, α, α_id, J, sub_J_ids, sup_α_ids) in bubble_functions
      passed == passed && w == w_prev+1
      w_prev = w

      passed = passed && all(α .≥ 0) && sum(α)==r-1 && length(α)==N &&
            α_id == bernstein_term_id(α) &&
            all( bernstein_term_id( [α[j]+Int(i==j) for j in eachindex(α)] ) == αpi_id
                for (i,αpi_id) in enumerate(sup_α_ids) )

      passed = passed && issorted(J) && length(J)==k+1 && (J ⊆ 1:N) &&
            all( _combination_index(J[J .≠ J[i]]) == Jsi_id for (i,Jsi_id) in enumerate(sub_J_ids) )
    end
    @test passed || (r, k, D, F, bubble_functions)
  end
  @test w_prev == binomial(r+k-1,k)*binomial(D+r,D-k)
end

for k in 0:D
  w_prev = 0
  for (F, bubble_functions) in PΛ_bubbles(r,k,D)
    @test issorted(F)  &&  (k ≤ length(F) ≤ N)  &&  (F ⊆ 1:N) || (k,D,F)

    passed = true
    for (w, α, α_id, J) in bubble_functions
      passed == passed && w == w_prev+1
      w_prev = w

      passed = passed && all(α .≥ 0) && sum(α)==r && length(α)==N &&
            α_id == bernstein_term_id(α)

      passed = passed && issorted(J) && length(J)==k && (J ⊆ 1:N)
    end
    @test passed || (r, k, D, F, bubble_functions)
  end
  @test w_prev == binomial(r+k,k)*binomial(D+r,D-k)
end


# TODO validation of precomputation in reference simplex

function _test_bases_in_ref_simplex(VD, T, r, k)
  b   = PmLambdaBasis(VD,T,r,k)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  bv  = PmLambdaBasis(VD,T,r,k,vertices)
  evaluate(bv,x)
  evaluate(Broadcasting(∇)(bv),x)
  evaluate(Broadcasting(∇∇)(bv),x)

  @test b.m == bv.m

  b   = PLambdaBasis(VD,T,r,k)
  evaluate(b,x)
  evaluate(Broadcasting(∇)(b),x)
  evaluate(Broadcasting(∇∇)(b),x)

  bv  = PLambdaBasis(VD,T,r,k,vertices)
  evaluate(bv,x)
  evaluate(Broadcasting(∇)(bv),x)
  evaluate(Broadcasting(∇∇)(bv),x)

  @test b.Ψ == bv.Ψ

  @test testvalue(typeof(b)) isa typeof(b)
end

T = Float64

# 0D                                           0D #
D = 0
vertices = (Point{D,T}(),)
x = [vertices[1]]
k = 0

_test_bases_in_ref_simplex(Val(D), T, r, k)

# 1D                                           1D #
D = 1
Pt = Point{D,T}
vertices = (Pt(0.),Pt(1.))
x = [xi for xi in vertices]

for k in 0:D
  _test_bases_in_ref_simplex(Val(D), T, r, k)
end


# 2D                                           2D #
D = 2
Pt = Point{D,T}
vertices = (Pt(0., 0),Pt(1.,0),Pt(0.,1.))
x = [xi for xi in vertices]

for k in 0:D
  _test_bases_in_ref_simplex(Val(D), T, r, k)
end

# 3D                                           3D #
D = 3
Pt = Point{D,T}
vertices = (Pt(0.,0,0),Pt(1.,0,0),Pt(0,1.,0),Pt(0,0,1.))
x = [xi for xi in vertices]

for k in 0:D
  _test_bases_in_ref_simplex(Val(D), T, r, k)
end

# 4D                                           4D #
D = 4
Pt = Point{D,T}
vertices = (Pt(0.,0,0,0),Pt(1.,0,0,0),Pt(0,1.,0,0),Pt(0,0,1.,0),Pt(0,0,0,1.))
x = [xi for xi in vertices]

for k in 0:D
  k == 2 && continue # no vector proxy
  _test_bases_in_ref_simplex(Val(D), T, r, k)
end

end # module
