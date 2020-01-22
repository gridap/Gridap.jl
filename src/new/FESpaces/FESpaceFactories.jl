
"""
"""
function FESpace(;kwargs...)

  constraint = _get_kwarg(:constraint,kwargs,nothing)
  reffe = _get_kwarg(:reffe,kwargs)
  @notimplementedif !isa(reffe,Symbol) "For the moment, reffe can only be a symbol"

  fespace = nothing

  if reffe in [:Lagrangian,:QLagrangian,:PLagrangian,:SLagrangian]

    fespace = _setup_lagrange_spaces(kwargs)

  else
    @notimplemented "For the moment only Lagrangian elements supported."
  end

  @assert fespace != nothing

  if constraint == nothing
    return fespace

  elseif constraint == :zeromean
    @notimplemented ":zeromean option not yet implemented"

  else
    @unreachable "Unknown constraint value $constraint"
  end

end

function _setup_lagrange_spaces(kwargs)

  reffe = _get_kwarg(:reffe,kwargs)
  model = _get_kwarg(:model,kwargs)
  labels = _get_kwarg(:labels,kwargs,nothing)
  conformity = _get_kwarg(:conformity,kwargs,true)
  diritags = _get_kwarg(:dirichlet_tags,kwargs,Int[])
  dirimasks = _get_kwarg(:dirichlet_masks,kwargs,nothing)
  order = _get_kwarg(:order,kwargs)

  polytopes = get_polytopes(model)
  trian = get_triangulation(model)

  T = _get_kwarg(:valuetype,kwargs,nothing)
  if T == nothing
    @unreachable "valuetype is a mandatory keyword argument in FESpace constructor for Lagrangian reference FEs"
  end

  if _is_reffe_lagrangian_compatible_with_polytopes(reffe,polytopes)

    if reffe == :SLagrangian
      _reffes = [SerendipityRefFE(T,p,order) for p in polytopes]
    else
      _reffes = [LagrangianRefFE(T,p,order) for p in polytopes]
    end

    if conformity in [false, :L2]

      s = "Strong dirichlet conditions cannot be imposed in discontinuous spaces for the moment"
      @notimplementedif diritags != Int[] s
      @notimplementedif dirimasks != nothing s

      return  DiscontinuousFESpace(_reffes,trian)

    elseif conformity in [true, :default, :H1]
      if labels == nothing
        return GradConformingFESpace(_reffes,model,diritags,dirimasks)
      else
        return GradConformingFESpace(_reffes,model,labels,diritags,dirimasks)
      end

    else
      s = "Conformity $conformity not implemented for $reffe reference FE on polytopes $(polytopes)"
      @unreachable s

    end

  else
    @notimplemented "PDisc spaces not yet implemented"
  end

end

function _get_kwarg(kwarg,kwargs)
  try
    return kwargs[kwarg]
  catch
    s = "The key-word argument $(kwarg) is mandatory in the FESpace constructor"
    error(s)
  end
end

function _get_kwarg(kwarg,kwargs,value)
  try
    return kwargs[kwarg]
  catch
    return value
  end
end

function _is_reffe_lagrangian_compatible_with_polytopes(reffe,polytopes)
  a = true
  for p in polytopes
    a = a && _is_reffe_lagrangian_compatible_with_polytope(reffe,p)
  end
  a
end

function _is_reffe_lagrangian_compatible_with_polytope(reffe,polytope)
  v = (reffe == :Lagrangian)
  v = v || (reffe == :QLagrangian && is_n_cube(polytope))
  v = v || (reffe == :SLagrangian && is_n_cube(polytope))
  v = v || (reffe == :PLagrangian && is_simplex(polytope))
  v
end
