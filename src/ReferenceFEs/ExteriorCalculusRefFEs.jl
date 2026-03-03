############################################
# (Proxied) Form valued FEEC reference FEs #
############################################

function ReferenceFE(p::Polytope{D},F::Symbol,r,k,T::Type;
  rotate_90=false, DG_calc=false, nodal::Bool=false, kwargs...) where D

  cart_prod = k in (0,D) && T <: MultiValue
  FEEC_space_definition_checks(Val(D),T,r,k,F,rotate_90,DG_calc; cart_prod)
  if !(k==D && F∈(:P,:S))
    @check r ≥ 1 "This exterior calculus FE starts at r=1, (F,r,k) = ($F,$r,$k)"
  end

  scalar = F==:S ? serendipity : lagrangian
  nodal && @check k in (0,D) "Nodal DOF only available for scalar FEs, i.e. k=0 or D, got k=$k, D=$D."
  scalar = nodal ? scalar : ModalScalar(scalar)

  # This logic should keep consistent with the Table in the ReferenceFEs documentation
  (name, order) = if F == :P⁻
    @check is_simplex(p)
    if     k == 0
      scalar, r
    elseif k == D
      scalar, r-1
    elseif k == 1 && !rotate_90
      nedelec1, r-1
    else # must be k = D-1 && rotate_90 = true if k = 1
      raviart_thomas, r-1
    end
  elseif F == :P
    @check is_simplex(p) #&& r>0
    if     k == 0
      #@check r>0 "P⁻Λᵏ starts at r=1"
      scalar, r
    elseif k == D
      scalar, r
    elseif k == 1 && !rotate_90
      nedelec2, r
    else # must be k = D-1 && rotate_90 = true if k = 1
      bdm, r
    end
  elseif F == :Q⁻
    @check is_n_cube(p) #&& r>0
    if     k == 0
      scalar, r
    elseif k == D
      scalar, r-1
    elseif k == 1 && !rotate_90
      nedelec1, r-1
    else # must be k = D-1 && rotate_90 = true if k = 1
      raviart_thomas, r-1
    end
  elseif F == :S
    @check is_n_cube(p) #&& r>0
    if     k == 0
      scalar, r
    elseif k == D
      if nodal
        return ReferenceFE(p,lagrangian,T,r; space=:P, kwargs...) # ℙᴰᵣ space ≡ SᵣΛᴰ
      else
        @notimplemented "Modal SrΛᴰ elements are not implemented, please use nodal ones (pass `nodal=true`)."
      end
    else
      @notimplemented "BDM on n-cubes and Serendipity Nedelec not implemented yet"
    #elseif k == 1 && !rotate_90
    #  nedelec2, r
    #else # must be k = D-1 && rotate_90 = true if k = 1
    #  bdm, r
    end
  else
    @unreachable
  end

  ReferenceFE(p,name,T,order; kwargs...)
end
