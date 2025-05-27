#############################################
# (Proxied) Form valued nD polynomial bases #
#############################################

"""
    FEECPolyBasis_trampoline(::Val{D},T,r,k,F::Symbol, PT=Monomial, vertices=nothing; kwargs...)

Return, if it is implemented, a polynomial basis for `D`-dimensional
(vector proxied) `k`-forms of polynomial order `r`. `T` is the scalar component
type, `PT<:Polynomial` the polynomial basis family, and `vertices` can sometimes
be given to build a basis in a physical polytope instead of the reference polytope.
"""
function FEECPolyBasis_trampoline(::Val{D},::Type{T},r,k,F::Symbol,::Type{PT}=Monomial,
    vertices=nothing; kwargs...) where {D,T,PT}

  @check F in (:P⁻,:P,:Q⁻,:S) "F must be either :P⁻,:P,:Q⁻ or :S"

  F == :P⁻ && PT == Bernstein && return PmLambdaBasis(Val(D),T,r,k,vertices; kwargs...)
  F == :P  && PT == Bernstein && return PLambdaBasis( Val(D),T,r,k,vertices; kwargs...)

  if k == 0
    # Scalar functions
    @notimplementedif r < 0
    if     F == :P⁻ || F == :P # Lagrange, ℙr space
      CartProdPolyBasis(PT,Val(D),T,r,_p_filter)
    elseif F == :Q⁻            # Lagrange, ℚr space
      CartProdPolyBasis(PT,Val(D),T,r,_q_filter)
    elseif F == :S             # Lagrange, 𝕊r space
      CartProdPolyBasis(PT,Val(D),T,r,_ser_filter)
    end


  elseif k == D
    # Scalar densities
    if     F == :P⁻ # Lagrange, ℙr₋1 space
      CartProdPolyBasis(PT,Val(D),T,r-1,_p_filter)
    elseif F == :P  # Lagrange, ℙr space
      CartProdPolyBasis(PT,Val(D),T,r,_p_filter)
    elseif F == :Q⁻ # Lagrange, ℚr₋1 space
      CartProdPolyBasis(PT,Val(D),T,r-1,_q_filter)
    elseif F == :S  # Serandipity Lagrange ≡ ℙr space
      CartProdPolyBasis(PT,Val(D),T,r,_p_filter)
    end


  elseif k == 1 # D > 1
    @notimplementedif D > 3
    V = VectorValue{D,T}

    if D == 2
      # TODO rot_90 for Curl/Div conform bases
      if     F == :P⁻ # Raviart-Thomas
        @notimplementedif PT ≠ Monomial
        PCurlGradBasis(PT,Val(D),T,r-1)
      elseif F == :P  # BDM on simplex
        CartProdPolyBasis(PT,Val(D),V,r,Polynomials._p_filter)
      elseif F == :Q⁻ # Raviart-Thomas
        QCurlGradBasis(PT,Val(D),T,r-1)
      elseif F == :S  # BDM on D-cubes
        @notimplemented
      end

    else # so D == 3
      if     F == :P⁻ # First kind Nedelec
        @notimplementedif PT ≠ Monomial
        PGradBasis(PT,Val(D),T,r-1)
      elseif F == :P  # Second kind Nedelec
        @notimplemented
      elseif F == :Q⁻ # First kind Nedelec
        QGradBasis(PT,Val(D),T,r-1)
      elseif F == :S  # "Serendipity second kind Nedelec" ?
        @notimplemented
      end
    end


  elseif k == 2 # D > 2
    @notimplementedif D > 3
    # D == 3
    if     F == :P⁻ # Raviart-Thomas
      PCurlGradBasis(PT,Val(D),T,r-1)
    elseif F == :P  # "3D BDM" ?
      @notimplemented
    elseif F == :Q⁻ # Raviart-Thomas
      QCurlGradBasis(PT,Val(D),T,r-1)
    elseif F == :S  # "3D Serandipity BDM" ?
      @notimplemented
    end

  else
    @unreachable
  end
end
