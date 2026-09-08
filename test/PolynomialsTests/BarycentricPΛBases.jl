module BarycentricPΛBasisTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Helpers
using ForwardDiff
using StaticArrays

using Gridap.Polynomials: _combination_index, bernstein_term_id

r = 3 # all possible bubble spaces are non empty

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


# Bases tests

using Gridap.Polynomials: _minusone_if_even_else_one, _findfirst_val_or_zero

function _test_testvalue(b, Bx, Gx, Hx)
  b0 = testvalue(b)
  @test b0 isa typeof(b)

  @test evaluate(b0,x)                   isa typeof(Bx)
  @test evaluate(Broadcasting(∇)(b0),x)  isa typeof(Gx)
  @test evaluate(Broadcasting(∇∇)(b0),x) isa typeof(Hx)
end

#############################################################################
# Former logic with analytical computation of the basis in the Reference FE #
#############################################################################

function _test_reference_basis(b,D,r,k)
  V = value_type(b)
  V <: Real && (V = VectorValue{1,V})

  if b isa BarycentricPmΛBasis
    LN = binomial(D+1,k)
    m = zero(MVector{LN,V})
    _compute_PmΛ_basis_reference_coefficients!(m,k,D,b._indices)
    @test all(@. norm(b.m - m) < 1.e-15)

  else       #BarycentricPΛBasis
    Ψ = zero(MVector{length(b),V})
    _compute_PΛ_basis_reference_form_coefficient!(Ψ,r,k,b._indices)
    @test all(@. norm(b.Ψ - Ψ) < 1.e-15)
  end
end

function _compute_PmΛ_basis_reference_coefficients!(m,k,D,indices)

  if iszero(k) # so V is scalar, no change of basis
    m .= 1
    return nothing
  end

  V = eltype(m)
  m_J = Mutable(V)(undef)
  @inbounds for (J_id, J) in enumerate(Polynomials._sorted_combinations(D+1,k))
    s = Int(isone(J[1]))
    for (I_id, I, I_sgn) in indices.components
      n = count(i-> (J[i]-1)∉I, (1+s):k)
      if iszero(n)
        p = _findfirst_val_or_zero(j-> (I[j]+1)∉J, 1, k)
        m_J[I_id] = I_sgn*_minusone_if_even_else_one(p+1)
      else
        m_J[I_id] = 0
      end
    end
    m[J_id] = m_J
  end
  nothing
end

function _compute_PΛ_basis_reference_form_coefficient!(Ψ,r,k,indices)

  """
      _hat_Ψ(r,::Val{k},α,F,I,J,T)::T

  BarycentricPΛBasis.Ψ matrix elements in the reference simplex, T is the scalar return type

  This is actually not faster than computing the matrices and the minors
  explicitely like when vertices are given, but might be usefull in case we want
  to compute these at compile time one day.
  """
  function _hat_Ψ(r,Vk::Val{k},α,F,I,J,::Type{T})::T where {T,k}
    @check sum(α) == r
    @check length(I) == length(J) == k

    iszero(k) && return one(T) # 0 forms

    _u(i::Int,F,I)   = Int(isone(F[1])) - Int(I[i]+1 in F)
    _u(F::Vector{Int},I,Vk) = ntuple(i->_u(i,F,I), Vk)
    _v(j::Int,α,J,r) = α[J[j]]/r
    _v(α::Vector{Int},J,r,Vk) = ntuple(j->_v(j,α,J,r), Vk)

    @inbounds begin

      s = Int(isone(J[1]))
      n = count(i-> (J[i]-1)∉I, (1+s):k)

      n > 1 && return 0. # rank M_IJ inferior to 2

      p = _findfirst_val_or_zero(j-> (I[j]+1)∉J, 1, k)

      if isone(n)        # rank M_IJ is 1
        m = _findfirst_val_or_zero(i-> (J[i]-1)∉I, (s+1), k)
        u_p, v_m = _u(p,F,I), _v(m,α,J,r)
        sgn = _minusone_if_even_else_one(m+p+1)
        iszero(s) && return sgn*u_p*v_m

        q = _findfirst_val_or_zero(j-> (I[j]+s)∉J, (p+1), k)
        u_q = _u(q,F,I)
        sgn *= _minusone_if_even_else_one(q+1)
        return sgn * v_m * (u_q - u_p)
      end

      u, v = _u(F,I,Vk), _v(α,J,r,Vk)
      if iszero(s)
        return 1 + sum( u .* v )
      else
        Ψ_IJ = one(T)
        sum_v = v[1]
        for l in 1:p-1
          vlp = v[l+1]
          sum_v += vlp
          Ψ_IJ += vlp*u[l]
        end
        for l in (p+1):k
          vl = v[l]
          sum_v += vl
          Ψ_IJ += vl*u[l]
        end
        sgn = _minusone_if_even_else_one(p+1)
        return sgn * (Ψ_IJ - u[p]*sum_v)
      end

    end
    @unreachable
  end

  Vk = Val(k)
  V = eltype(Ψ)
  T = eltype(V)
  Ψw = Mutable(V)(undef)

  iszero(r) && return _order_0_Ψ!(Ψ)

  @inbounds for (F, bubble_functions) in indices.bubbles
    for (w, α, _, J) in bubble_functions
      for (I_id, I, I_sgn) in indices.components
        Ψw[I_id] = I_sgn * _hat_Ψ(r,Vk,α,F,I,J,T)
      end
      Ψ[w] = Ψw
    end
  end
  nothing
end

##############################
# testing function and tests #
##############################

function _test_basis(VD::Val{D}, T, r, k, vertices) where D
  for PΛB in (BarycentricPmΛBasis, BarycentricPΛBasis)
    b   = PΛB(VD,T,r,k)
    @test contains(sprint(show, MIME"text/plain"(), b._indices), "PᵣΛᵏ(△ᴰ) basis indices, r=$r k=$k D=$D")
    @test_nowarn print_indices(b,IOBuffer())
    @test get_orders(b) == tfill(r,Val(D))

    b2  = PΛB(VD,T,r,k; indices=b._indices) # indices recycling
    @test b == b2
    @test b2._indices == b._indices

    faces = [bubble[1] for bubble in get_bubbles(b)] # bubble space selection
    b2  = PΛB(b, faces...)
    @test b == b2

    _test_reference_basis(b,D,r,k)

    Bx = evaluate(b,x)
    Gx = evaluate(Broadcasting(∇)(b),x)
    Hx = evaluate(Broadcasting(∇∇)(b),x)
    _test_testvalue(b, Bx, Gx, Hx)

    bv  = PΛB(VD,T,r,k,vertices)
    Bx = evaluate(bv,x)
    Gx = evaluate(Broadcasting(∇)(bv),x)
    Hx = evaluate(Broadcasting(∇∇)(bv),x)
    _test_testvalue(bv, Bx, Gx, Hx)
  end
end

T = Float64

# 0D                                           0D #
D = 0
vertices = [Point{D,T}()]
x = [vertices[1]]
k = 0

_test_basis(Val(D), T, r, k, vertices)

# 1D                                           1D #
D = 1
Pt = Point{D,T}
vertices = [Pt(0.),Pt(1.)]
x = [xi for xi in vertices]

for k in 0:D
  _test_basis(Val(D), T, r, k, vertices)
end


# 2D                                           2D #
D = 2
Pt = Point{D,T}
vertices = [Pt(0., 0),Pt(1.,0),Pt(0.,1.)]
x = [xi for xi in vertices]

for k in 0:D
  _test_basis(Val(D), T, r, k, vertices)
end

# 3D                                           3D #
D = 3
Pt = Point{D,T}
vertices = [Pt(0.,0,0),Pt(1.,0,0),Pt(0,1.,0),Pt(0,0,1.)]
x = [xi for xi in vertices]

for k in 0:D
  _test_basis(Val(D), T, r, k, vertices)
end

# 4D                                           4D #
D = 4
Pt = Point{D,T}
vertices = [Pt(0.,0,0,0),Pt(1.,0,0,0),Pt(0,1.,0,0),Pt(0,0,1.,0),Pt(0,0,0,1.)]
x = [xi for xi in vertices]

for k in 0:D
  k == 2 && continue # no vector proxy
  _test_basis(Val(D), T, r, k, vertices)
end

end # module
