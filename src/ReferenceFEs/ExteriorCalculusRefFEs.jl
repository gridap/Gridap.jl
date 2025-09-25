############################################
# (Proxied) Form valued FEEC reference FEs #
############################################

function ReferenceFE(p::Polytope{D},F::Symbol,r,k,T::Type;
  rotate_90=false, DG_calc=false, kwargs...) where D

  FEEC_space_definition_checks(Val(D),T,r,k,F,rotate_90,DG_calc)
  if !(k==D && F∈(:P,:S))
    @check r ≥ 1 "This exterior calculus FE starts at r=1, (F,r,k) = ($F,$r,$k)"
  end

  # This logic should keep consistent with the Table in the ReferenceFEs documentation
  (name, order) = if F == :P⁻
    @check is_simplex(p)
    if     k == 0
      Lagrangian(), r
    elseif k == D
      Lagrangian(), r-1
    elseif k == 1 && !rotate_90
      Nedelec{1}(), r-1
    else # must be k = D-1 && rotate_90 = true if k = 1
      RaviartThomas(), r-1
    end
  elseif F == :P
    @check is_simplex(p) #&& r>0
    if     k == 0
      #@check r>0 "P⁻Λᵏ starts at r=1"
      Lagrangian(), r
    elseif k == D
      Lagrangian(), r
    elseif k == 1 && !rotate_90
      Nedelec{2}(), r
    else # must be k = D-1 && rotate_90 = true if k = 1
      BDM(), r
    end
  elseif F == :Q⁻
    @check is_n_cube(p) #&& r>0
    if     k == 0
      Lagrangian(), r
    elseif k == D
      Lagrangian(), r-1
    elseif k == 1 && !rotate_90
      Nedelec{1}(), r-1
    else # must be k = D-1 && rotate_90 = true if k = 1
      RaviartThomas(), r-1
    end
  elseif F == :S
    @check is_n_cube(p) #&& r>0
    if     k == 0
      Serendipity(), r
    elseif k == D
      return ReferenceFE(p,Lagrangian(),T,r; space=:P) # ℙᴰᵣ space ≡ SᵣΛᴰ
    else
      @notimplemented "BDM on n-cubes and Serendipity Nedelec not implemented yet"
    #elseif k == 1 && !rotate_90
    #  Nedelec2(), r
    #else # must be k = D-1 && rotate_90 = true if k = 1
    #  BDM(), r
    end
  else
    @unreachable
  end

  ReferenceFE(p,name,order; kwargs...)
end
