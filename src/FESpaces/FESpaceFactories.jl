
"""
"""
function FESpace(;kwargs...)

  constraint = _get_kwarg(:constraint,kwargs,nothing)
  reffe = _get_kwarg(:reffe,kwargs)
  @notimplementedif !isa(reffe,Symbol) "For the moment, reffe can only be a symbol"

  fespace = nothing

  if reffe in [:Lagrangian,:QLagrangian,:PLagrangian,:SLagrangian]

    fespace = _setup_lagrange_spaces(kwargs)

  elseif reffe == :RaviartThomas

    fespace = _setup_hdiv_space(kwargs)

  elseif reffe == :Nedelec

    fespace = _setup_hcurl_space(kwargs)

  else
    @notimplemented "Unsupported reffe $reffe."
  end

  @assert fespace != nothing

  if constraint == nothing
    return fespace

  elseif constraint == :zeromean
    model = _get_kwarg(:model,kwargs)
    order = _get_kwarg(:order,kwargs)
    trian = get_triangulation(model)
    quad = CellQuadrature(trian,order)
    return ZeroMeanFESpace(fespace,trian,quad)

  else
    @unreachable "Unknown constraint value $constraint"
  end

end

"""
"""
function TestFESpace(;kwargs...)
  FESpace(;kwargs...)
end

function _setup_hdiv_space(kwargs)

  reffe = _get_kwarg(:reffe,kwargs)
  model = _get_kwarg(:model,kwargs)
  labels = _get_kwarg(:labels,kwargs,get_face_labeling(model))
  conformity = _get_kwarg(:conformity,kwargs,true)
  diritags = _get_kwarg(:dirichlet_tags,kwargs,Int[])
  order = _get_kwarg(:order,kwargs,nothing)
  Tf = _get_kwarg(:valuetype,kwargs,VectorValue{1,Float64})
  T = eltype(Tf)

  if order == nothing
    @unreachable "order is a mandatory keyword argument in FESpace constructor for RaviartThomas reference FEs"
  end

  polytopes = get_polytopes(model)
  reffes = [RaviartThomasRefFE(T,p,order) for p in polytopes]

  if conformity in [true, :default, :HDiv, :Hdiv]
      V =  DivConformingFESpace(reffes,model,labels,diritags)
  else
    s = "Conformity $conformity not implemented for $reffe reference FE on polytopes $(polytopes...)"
    @unreachable s
  end

  V
end

function _setup_hcurl_space(kwargs)

  reffe = _get_kwarg(:reffe,kwargs)
  model = _get_kwarg(:model,kwargs)
  labels = _get_kwarg(:labels,kwargs,get_face_labeling(model))
  conformity = _get_kwarg(:conformity,kwargs,true)
  diritags = _get_kwarg(:dirichlet_tags,kwargs,Int[])
  order = _get_kwarg(:order,kwargs,nothing)
  Tf = _get_kwarg(:valuetype,kwargs,VectorValue{1,Float64})
  T = eltype(Tf)

  if order == nothing
    @unreachable "order is a mandatory keyword argument in FESpace constructor for Nedelec reference FEs"
  end

  polytopes = get_polytopes(model)
  reffes = [NedelecRefFE(T,p,order) for p in polytopes]

  if conformity in [true, :default, :HCurl, :Hcurl]
      V =  CurlConformingFESpace(reffes,model,labels,diritags)
  else
    s = "Conformity $conformity not implemented for $reffe reference FE on polytopes $(polytopes...)"
    @unreachable s
  end

  V
end

function _setup_lagrange_spaces(kwargs)

  conformity = _get_kwarg(:conformity,kwargs,true)
  reffe = _get_kwarg(:reffe,kwargs)
  order = _get_kwarg(:order,kwargs)
  T = _get_kwarg(:valuetype,kwargs,nothing)
  diritags = _get_kwarg(:dirichlet_tags,kwargs,Int[])
  dirimasks = _get_kwarg(:dirichlet_masks,kwargs,nothing)
  labels = _get_kwarg(:labels,kwargs,nothing)
  model = _get_kwarg(:model,kwargs,nothing)

  if T == nothing
    @unreachable "valuetype is a mandatory keyword argument in FESpace constructor for Lagrangian reference FEs"
  end

  if conformity in [false, :L2]

    s = "Strong dirichlet conditions cannot be imposed in discontinuous spaces for the moment"
    @notimplementedif diritags != Int[] s
    @notimplementedif dirimasks != nothing s

    _trian = _get_kwarg(:triangulation,kwargs,nothing)
    if _trian == nothing
      if model == nothing
        @unreachable "either a model or a triangulation has to be provided for building a discontinuous Lagrangian space"
      end
      trian = get_triangulation(model)
    else
      if model != nothing
        @unreachable "either a model or a triangulation BUT NOT BOTH has to be provided for building a discontinuous Lagrangian space"
      end
      trian = _trian
    end

    polytopes = [get_polytope(r) for r in get_reffes(trian)]

    if _is_reffe_lagrangian_compatible_with_polytopes(reffe,polytopes)
      if reffe == :SLagrangian
        _reffes = [SerendipityRefFE(T,p,order) for p in polytopes]
      else
        _reffes = [LagrangianRefFE(T,p,order) for p in polytopes]
      end
    else
      if reffe == :PLagrangian
        _reffes = [PDiscRefFE(T,p,order) for p in polytopes]
      else
        @unreachable "Not possible to use a $reffe reffe on polytopoes $(polytopes...)"
      end
    end

    return  DiscontinuousFESpace(_reffes,trian)

  elseif conformity in [true, :default, :H1, :C0]

    if model == nothing
      @unreachable "model is a mandatory keyword argument in FESpace constructor for conforming Lagrangian reference FEs"
    end
    polytopes = get_polytopes(model)
    trian = get_triangulation(model)
    if ! _is_reffe_lagrangian_compatible_with_polytopes(reffe,polytopes)
      s = "Conformity $conformity not possible for $reffe reference FE on polytopes $(polytopes...)"
      @unreachable s
    end
    if reffe == :SLagrangian
      _reffes = [SerendipityRefFE(T,p,order) for p in polytopes]
    else
      _reffes = [LagrangianRefFE(T,p,order) for p in polytopes]
    end
    if labels == nothing
      return GradConformingFESpace(_reffes,model,diritags,dirimasks)
    else
      return GradConformingFESpace(_reffes,model,labels,diritags,dirimasks)
    end

  else
    s = "Conformity $conformity not implemented for lagrangian reference FEs"
    @unreachable s
  end

end

#function _setup_lagrange_spaces(kwargs)
#
#  reffe = _get_kwarg(:reffe,kwargs)
#  model = _get_kwarg(:model,kwargs)
#  labels = _get_kwarg(:labels,kwargs,nothing)
#  conformity = _get_kwarg(:conformity,kwargs,true)
#  diritags = _get_kwarg(:dirichlet_tags,kwargs,Int[])
#  dirimasks = _get_kwarg(:dirichlet_masks,kwargs,nothing)
#  order = _get_kwarg(:order,kwargs)
#
#  polytopes = get_polytopes(model)
#  trian = get_triangulation(model)
#
#  T = _get_kwarg(:valuetype,kwargs,nothing)
#  if T == nothing
#    @unreachable "valuetype is a mandatory keyword argument in FESpace constructor for Lagrangian reference FEs"
#  end
#
#  if _is_reffe_lagrangian_compatible_with_polytopes(reffe,polytopes)
#
#    if reffe == :SLagrangian
#      _reffes = [SerendipityRefFE(T,p,order) for p in polytopes]
#    else
#      _reffes = [LagrangianRefFE(T,p,order) for p in polytopes]
#    end
#
#    if conformity in [false, :L2]
#
#      s = "Strong dirichlet conditions cannot be imposed in discontinuous spaces for the moment"
#      @notimplementedif diritags != Int[] s
#      @notimplementedif dirimasks != nothing s
#
#      return  DiscontinuousFESpace(_reffes,trian)
#
#    elseif conformity in [true, :default, :H1, :C0]
#      if labels == nothing
#        return GradConformingFESpace(_reffes,model,diritags,dirimasks)
#      else
#        return GradConformingFESpace(_reffes,model,labels,diritags,dirimasks)
#      end
#
#    else
#      s = "Conformity $conformity not implemented for $reffe reference FE on polytopes $(polytopes...)"
#      @unreachable s
#
#    end
#
#  elseif reffe == :PLagrangian
#
#      if conformity in [false, :L2]
#
#        _reffes = [PDiscRefFE(T,p,order) for p in polytopes]
#        return  DiscontinuousFESpace(_reffes,trian)
#
#      else
#
#        @unreachable "Conformity $conformity not possible for $reffe reference FE on $(polytopes...)"
#
#      end
#
#  else
#
#    @notimplemented "Reference element $reffe not implemented on $(polytopes...)"
#
#  end
#
#end

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
