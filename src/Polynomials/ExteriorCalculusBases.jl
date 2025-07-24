#############################################
# (Proxied) Form valued nD polynomial bases #
#############################################

_ensure_hierarchical(PT) = !isHierarchical(PT) && @unreachable "Polynomial family must be hierarchical, got $PT."

function _default_poly_type(F)
  F ∈ (:P, :P⁻) && return Bernstein
  F ∈ (:Q⁻,:S)  && return Legendre
  Monomial
end

"""
    FEEC_poly_basis(::Val{D},T,r,k,F::Symbol, PT=_default_poly_type(F); kwargs...)

"Factory for polynomial basis of Finite Element Exterior Calculus spaces"

Return, if it is implemented, a polynomial basis for the space  `FᵣΛᵏ` in
dimension `D`, with `T` the scalar component type and `PT<:Polynomial` the
polynomial basis family.

The default `PT` is `Bernstein` on simplices and `Legendre` on D-cubes.

# Arguments
- `D`: spatial dimension
- `T::Type`: scalar components type
- `r::Int`: polynomial order
- `k::Int`: form order
- `F::Symbol`: family, i.e. `:P⁻`, `:P`, `:Q⁻` or `:S`
### kwargs
- `rotate_90::Bool`: only if `D`=2 and `k`=1, tells to use the vector proxy corresponding to div conform function instead of curl conform ones.
- `vertices=nothing`: for `PT=Bernstein` bases on simplices (`F = :P` or `:P⁻`), the basis is defined on the simplex defined by `vertices` instead of the reference one.
""" # document diff_geo_calculus_style once its implemented
function FEEC_poly_basis(::Val{D},::Type{T},r,k,F::Symbol,PT=_default_poly_type(F);
    diff_geo_calculus_style=false, rotate_90=false, vertices=nothing) where {D,T}

  @assert PT <: Polynomial

  # these call FEEC_space_definition_checks internally
  F == :P⁻ && PT == Bernstein && return PmLambdaBasis(Val(D),T,r,k,vertices)
  F == :P  && PT == Bernstein && return PLambdaBasis( Val(D),T,r,k,vertices)

  FEEC_space_definition_checks(Val(D), T, r, k, F, rotate_90, diff_geo_calculus_style)
  @notimplementedif diff_geo_calculus_style # This ensures 0≤k≤D≤3

  if k == 0
    # Scalar H1 conforming functions
    @notimplementedif r < 0
    if     F == :P⁻ || F == :P # Lagrange, ℙr space
      _ensure_hierarchical(PT)
      CartProdPolyBasis(PT,Val(D),T,r,_p_filter)
    elseif F == :Q⁻            # Lagrange, ℚr space
      CartProdPolyBasis(PT,Val(D),T,r,_q_filter)
    elseif F == :S             # Lagrange, 𝕊r space
      _ensure_hierarchical(PT)
      CartProdPolyBasis(PT,Val(D),T,r,_ser_filter)
    end


  elseif k == D
    # Scalar L2 conforming densities
    if     F == :P⁻ # Lagrange, ℙr₋1 space
      _ensure_hierarchical(PT)
      CartProdPolyBasis(PT,Val(D),T,r-1,_p_filter)
    elseif F == :P  # Lagrange, ℙr space
      _ensure_hierarchical(PT)
      CartProdPolyBasis(PT,Val(D),T,r,_p_filter)
    elseif F == :Q⁻ # Lagrange, ℚr₋1 space
      CartProdPolyBasis(PT,Val(D),T,r-1,_q_filter)
    elseif F == :S  # Serendipity Lagrange ≡ ℙr space
      _ensure_hierarchical(PT)
      CartProdPolyBasis(PT,Val(D),T,r,_p_filter)
    end


  elseif k == 1 # and D > 1
    V = VectorValue{D,T}
    if D == 2
      if     F == :P⁻
        if rotate_90 # Raviart-Thomas
          # former PCurlGradBasis(PT,Val(D),T,r-1)
          RaviartThomasPolyBasis{D}(PT, T, r-1)
        else         # Nedelec
          # former PGradBasis(PT,Val(D),T,r-1)
          _ensure_hierarchical(PT)
          NedelecPolyBasisOnSimplex{D}(PT,T,r-1)
        end
      elseif F == :P  # BDM on simplex
        CartProdPolyBasis(PT,Val(D),V,r,Polynomials._p_filter) # rotation not needed
      elseif F == :Q⁻
        if rotate_90 # Raviart-Thomas
          # former QCurlGradBasis(PT,Val(D),T,r-1)
          m = Tuple( r-1 + (i==j ? 1 : 0) for i in 1:D, j in 1:D )
          orders = SMatrix{D,D,Int}(m)
          CompWiseTensorPolyBasis{D}(PT, V, orders)
        else         # Nedelec
          # former QGradBasis(PT,Val(D),T,r-1)
          m = Tuple( r-1 + (i==j ? 0 : 1) for i in 1:D, j in 1:D )
          orders = SMatrix{D,D,Int}(m)
          CompWiseTensorPolyBasis{D}(PT, V, orders)
        end
      elseif F == :S  # BDM on D-cubes
        @notimplemented
      end

    elseif D == 3
      if     F == :P⁻ # First kind Nedelec
        # former PGradBasis(PT,Val(D),T,r-1)
        _ensure_hierarchical(PT)
        NedelecPolyBasisOnSimplex{D}(PT,T,r-1)
      elseif F == :P  # Second kind Nedelec
        #@notimplemented
        _ensure_hierarchical(PT)
        CartProdPolyBasis(PT,Val(D),V,r,Polynomials._p_filter) # rotation not needed
      elseif F == :Q⁻ # First kind Nedelec
        # former QGradBasis(PT,Val(D),T,r-1)
        m = Tuple( r-1 + (i==j ? 0 : 1) for i in 1:D, j in 1:D )
        orders = SMatrix{D,D,Int}(m)
        CompWiseTensorPolyBasis{D}(PT, V, orders)
      elseif F == :S  # "Serendipity second kind Nedelec" ?
        @notimplemented
      end

    else # D > 3
      @notimplemented
    end


  elseif k == 2 # D > 2
    @notimplementedif D > 3
    V = VectorValue{D,T}
    # D == 3
    if     F == :P⁻ # Raviart-Thomas
      # former PCurlGradBasis(PT,Val(D),T,r-1)
      RaviartThomasPolyBasis{D}(PT, T, r-1)
    elseif F == :P  # "3D BDM" ?
      _ensure_hierarchical(PT)
      CartProdPolyBasis(PT,Val(D),V,r,Polynomials._p_filter) # rotation not needed
    elseif F == :Q⁻ # Raviart-Thomas
      # former QCurlGradBasis(PT,Val(D),T,r-1)
      m = Tuple( r-1 + (i==j ? 1 : 0) for i in 1:D, j in 1:D )
      orders = SMatrix{D,D,Int}(m)
      CompWiseTensorPolyBasis{D}(PT, V, orders)
    elseif F == :S  # "3D Serendipity BDM" ?
      @notimplemented
    end

  else
    @unreachable
  end
end
