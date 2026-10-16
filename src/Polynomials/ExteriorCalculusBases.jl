#############################################
# (Proxied) Form valued nD polynomial bases #
#############################################

_ensure_hierarchical(PT) = !isHierarchical(PT) && @unreachable "Polynomial family must be hierarchical, got $PT."

function _default_poly_type(F)
  F ∈ (:P, :P⁻) && return Bernstein
  F ∈ (:Q⁻,:S)  && return ModalC0
  Monomial
end

"""
    FEEC_poly_basis(::Val{D},T,r,k, F::Symbol, PT=_default_poly_type(F); kwargs...)

"Factory for polynomial basis of Finite Element Exterior Calculus spaces"

Return, if it is implemented, a polynomial basis for the space  ``\\mathrm{F}ᵣΛᵏ`` in
dimension `D`, with  `PT<:Polynomial` the polynomial basis family.

`T` may be a `ExteriorFormValue{k,D}`, in which case the basis is `T` valued.
Else, `T` is the scalar type of the vector proxy of the space, except if
`cart_prod=true` is given, see kwargs below.

The default `PT` is `Bernstein` on simplices and `ModalC0` on D-cubes.

# Arguments
- `D`: spatial dimension
- `T::Type`: an `ExteriorFormValue` or the scalar component type, except if `cart_prod=true`, see kwargs
- `r::Int`: polynomial order
- `k::Int`: form order
- `F::Symbol`: family, i.e. `:P⁻`, `:P`, `:Q⁻` or `:S`
### kwargs
- `rotate_90::Bool`: only if `D`=2 and `k`=1, tells to use the vector proxy corresponding to div conform function instead of curl conform ones.
- `vertices=nothing`: for `PT=Bernstein` bases on simplices (`F = :P` or `:P⁻`), the basis is defined on the simplex defined by `vertices` instead of the reference one.
- `cart_prod=false`: for `k`=0 or `k=D` and if `T` is not a `ExteriorFormValue`, `true` authorise `T` to be a tensor type for Cartesian product space of the scalar space.
"""
function FEEC_poly_basis(::Val{D},::Type{T},r,k,F::Symbol,PT=_default_poly_type(F);
    rotate_90=false, vertices=nothing, cart_prod=false) where {D,T}

  @assert PT <: Polynomial

  # these call FEEC_space_definition_checks internally
  if !cart_prod
    F == :P⁻ && PT == Bernstein && return BarycentricPmΛBasis(Val(D),T,r,k,vertices; rotate_90)
    F == :P  && PT == Bernstein && return BarycentricPΛBasis( Val(D),T,r,k,vertices; rotate_90)
  end


  FEEC_space_definition_checks(Val(D), T, r, k, F, rotate_90; cart_prod)
  # This also ensures 0≤k≤D≤3 on the proxied path
  DG_calc = T <: ExteriorFormValue
  @notimplementedif DG_calc && !(k in (0,D))
  V = make_concretetype(T)

  if k == 0
    # Scalar H1 conforming functions
    @notimplementedif r < 0
    if     F == :P⁻ || F == :P # Lagrange, 𝓟r space
      PT == Bernstein && return BernsteinBasisOnSimplex{D}(V,r,vertices) # only if cart_prod or DG_calc
      _ensure_hierarchical(PT)
      CartProdPolyBasis(PT,Val(D),V,r,_p_filter)
    elseif F == :Q⁻            # Lagrange, 𝓠r space
      CartProdPolyBasis(PT,Val(D),V,r,_q_filter)
    elseif F == :S             # Lagrange, 𝓢r space
      PT==ModalC0 || _ensure_hierarchical(PT)
      CartProdPolyBasis(PT,Val(D),V,r,_ser_filter)
    end


  elseif k == D
    # Scalar L2 conforming densities
    if     F == :P⁻ # Lagrange, 𝓟r₋1 space
      PT == Bernstein && return BernsteinBasisOnSimplex{D}(V,r-1,vertices) # only if cart_prod or DG_calc
      _ensure_hierarchical(PT)
      CartProdPolyBasis(PT,Val(D),V,r-1,_p_filter)
    elseif F == :P  # Lagrange, 𝓟r space
      PT == Bernstein && return BernsteinBasisOnSimplex{D}(V,r,vertices) # only if cart_prod or DG_calc
      _ensure_hierarchical(PT)
      CartProdPolyBasis(PT,Val(D),V,r,_p_filter)
    elseif F == :Q⁻ # Lagrange, 𝓠r₋1 space
      CartProdPolyBasis(PT,Val(D),V,r-1,_q_filter)
    elseif F == :S  # Serendipity Lagrange ≡ 𝓟r space
      PT == Bernstein && return BernsteinBasisOnSimplex{D}(V,r,vertices) # only if cart_prod or DG_calc
      _ensure_hierarchical(PT)
      CartProdPolyBasis(PT,Val(D),V,r,_p_filter)
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
          NedelecPolyBasisOnSimplex{D}(PT,T,r-1)
        end
      elseif F == :P  # BDM on simplex, cannot be Bernstein
        _ensure_hierarchical(PT)
        CartProdPolyBasis(PT,Val(D),V,r,Polynomials._p_filter) # rotation not needed
      elseif F == :Q⁻
        if rotate_90 # Raviart-Thomas
          # former QCurlGradBasis(PT,Val(D),T,r-1)
          orders = [ r-1 + (i==j ? 1 : 0) for i in 1:D, j in 1:D ]
          CompWiseTensorPolyBasis{D}(PT, V, orders)
        else         # Nedelec
          # former QGradBasis(PT,Val(D),T,r-1)
          orders = [ r-1 + (i==j ? 0 : 1) for i in 1:D, j in 1:D ]
          CompWiseTensorPolyBasis{D}(PT, V, orders)
        end
      elseif F == :S  # BDM on D-cubes
        @notimplemented
      end

    elseif D == 3
      if     F == :P⁻ # First kind Nedelec, cannot be Bernstein
        # former PGradBasis(PT,Val(D),T,r-1)
        NedelecPolyBasisOnSimplex{D}(PT,T,r-1)
      elseif F == :P  # Second kind Nedelec, cannot be Bernstein
        _ensure_hierarchical(PT)
        CartProdPolyBasis(PT,Val(D),V,r,Polynomials._p_filter) # rotation not needed
      elseif F == :Q⁻ # First kind Nedelec
        # former QGradBasis(PT,Val(D),T,r-1)
        orders = [ r-1 + (i==j ? 0 : 1) for i in 1:D, j in 1:D ]
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
    if     F == :P⁻ # Raviart-Thomas, cannot be Bernstein
      # former PCurlGradBasis(PT,Val(D),T,r-1)
      RaviartThomasPolyBasis{D}(PT, T, r-1)
    elseif F == :P  # "3D BDM" ?  cannot be Bernstein
      _ensure_hierarchical(PT)
      CartProdPolyBasis(PT,Val(D),V,r,Polynomials._p_filter) # rotation not needed
    elseif F == :Q⁻ # Raviart-Thomas
      # former QCurlGradBasis(PT,Val(D),T,r-1)
      orders = [ r-1 + (i==j ? 1 : 0) for i in 1:D, j in 1:D ]
      CompWiseTensorPolyBasis{D}(PT, V, orders)
    elseif F == :S  # "3D Serendipity BDM" ?
      @notimplemented
    end

  else
    @unreachable
  end
end
