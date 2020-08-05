
"""
"""
function FESpace(;kwargs...)

  constraint = _get_kwarg(:constraint,kwargs,nothing)
  reffe = _get_kwarg(:reffe,kwargs)

  fespace = nothing

  if isa(reffe,ReferenceFE)

    fespace = _setup_generic_space(kwargs)

  elseif reffe in [:Lagrangian,:QLagrangian,:PLagrangian,:SLagrangian]

    fespace = _setup_lagrange_spaces(kwargs)

  elseif reffe == :RaviartThomas

    fespace = _setup_hdiv_space(kwargs)

  elseif reffe == :Nedelec

    fespace = _setup_hcurl_space(kwargs)

  else
    @notimplemented "Unsupported reffe $reffe."
  end

  @assert fespace != nothing

  restricted_at = _get_restricted_triangulation(kwargs)
  if restricted_at == nothing
    _fespace = fespace
  else
    @assert isa(restricted_at,RestrictedTriangulation)
    _fespace = ExtendedFESpace(fespace,restricted_at)
  end

  if constraint == nothing
    return _fespace

  elseif constraint == :zeromean
    order = _get_kwarg(:order,kwargs)
    model = _get_kwarg(:model,kwargs,nothing)
    zeromean_trian = _get_kwarg(:zeromean_trian,kwargs,nothing)
    zeromean_quad = _get_kwarg(:zeromean_quad,kwargs,nothing)
    if zeromean_trian !== nothing
      trian = zeromean_trian
    elseif model === nothing
      trian = _get_kwarg(:triangulation,kwargs)
    else
      trian = get_triangulation(model)
    end
    if zeromean_quad === nothing
      quad = CellQuadrature(trian,order)
    else
      quad = zeromean_quad
    end
    return ZeroMeanFESpace(_fespace,trian,quad)

  else
    @unreachable "Unknown constraint value $constraint"
  end


end

function _get_restricted_triangulation(kwargs)
  model  = _get_kwarg(:model,kwargs,nothing)
  _trian = _get_kwarg(:triangulation,kwargs,nothing)
  if _trian == nothing
    if model == nothing
      @unreachable "either a model or a triangulation has to be provided for building a FESpace"
    end
    trian = get_triangulation(model)
  else
    if model != nothing
      @unreachable "either a model or a triangulation BUT NOT BOTH has to be provided for building a FESpace"
    end
    trian = _trian
  end
  _get_restricted_triangulation_or_nothing(trian)
end

function  _get_restricted_triangulation_or_nothing(trian)
  nothing
end

function  _get_restricted_triangulation_or_nothing(trian::RestrictedTriangulation)
  trian
end

"""
"""
function TestFESpace(;kwargs...)
  FESpace(;kwargs...)
end

function _setup_generic_space(kwargs)

  reffe = _get_kwarg(:reffe,kwargs)
  model = _get_kwarg(:model,kwargs)
  labels = _get_kwarg(:labels,kwargs,get_face_labeling(model))
  conformity = _get_kwarg(:conformity,kwargs,true)
  diritags = _get_kwarg(:dirichlet_tags,kwargs,Int[])
  order = _get_kwarg(:order,kwargs,nothing)
  dofspace = _get_kwarg(:dof_space,kwargs,:reference)
  ( dofspace == :reference ? true : false )

  is_ref = (dofspace==:reference)

  Tf = _get_kwarg(:valuetype,kwargs,VectorValue{1,Float64})
  T = eltype(Tf)


  polytopes = get_polytopes(model)
  if length(polytopes) > 2 || polytopes[1] != get_polytope(reffe)
    @unreachable "RefFE and geometrical model are not compatible in generic FESpace constructor"
  end

  reffes = [ reffe ]

  if conformity in [true, :default]
    conf = get_default_conformity(reffe)
    V =  ConformingFESpace(reffes,conf,model,labels,diritags,nothing,is_ref)
  elseif isa(conformity,Conformity)
    V =  ConformingFESpace(reffes,conformity,model,labels,diritags,nothing,is_ref)
  else
    s = "Unsuported conformity = $conformity"
    @unreachable s
  end

  V
end

function _setup_hdiv_space(kwargs)

  reffe = _get_kwarg(:reffe,kwargs)
  model = _get_kwarg(:model,kwargs)
  labels = _get_kwarg(:labels,kwargs,get_face_labeling(model))
  conformity = _get_kwarg(:conformity,kwargs,true)
  diritags = _get_kwarg(:dirichlet_tags,kwargs,Int[])
  order = _get_kwarg(:order,kwargs,nothing)
  dofspace = _get_kwarg(:dof_space,kwargs,:reference)
  ( dofspace == :reference ? true : false )

  is_ref = (dofspace==:reference)

  Tf = _get_kwarg(:valuetype,kwargs,VectorValue{1,Float64})
  T = eltype(Tf)

  if order == nothing
    @unreachable "order is a mandatory keyword argument in FESpace constructor for RaviartThomas reference FEs"
  end

  polytopes = get_polytopes(model)
  reffes = [RaviartThomasRefFE(T,p,order) for p in polytopes]

  if conformity in [true, :default, :HDiv, :Hdiv]
    V =  ConformingFESpace(reffes,DivConformity(),model,labels,diritags,nothing,is_ref)
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
  dofspace = _get_kwarg(:dof_space,kwargs,:reference)
  Tf = _get_kwarg(:valuetype,kwargs,VectorValue{1,Float64})
  T = eltype(Tf)

  is_ref = dofspace==:reference

  if order == nothing
    @unreachable "order is a mandatory keyword argument in FESpace constructor for Nedelec reference FEs"
  end

  polytopes = get_polytopes(model)
  reffes = [NedelecRefFE(T,p,order) for p in polytopes]

  if conformity in [true, :default, :HCurl, :Hcurl]
    V =  ConformingFESpace(reffes,CurlConformity(),model,labels,diritags,nothing,is_ref)
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
  # dirimasks = nothing
  dofspace = _get_kwarg(:dof_space,kwargs,:reference)
  labels = _get_kwarg(:labels,kwargs,nothing)
  # labels = nothing
  model = _get_kwarg(:model,kwargs,nothing)

  is_ref = (dofspace==:reference)

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
        _reffes = [LagrangianRefFE(T,p,order;space=:P) for p in polytopes]
      else
        @unreachable "Not possible to use a $reffe reffe on polytopes $(polytopes...)"
      end
    end

    return  DiscontinuousFESpace(_reffes,trian,is_ref)

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
      return GradConformingFESpace(_reffes,model,diritags,dirimasks,is_ref)
    else
      return ConformingFESpace(_reffes,GradConformity(),model,labels,diritags,dirimasks,is_ref)
    end

  else
    s = "Conformity $conformity not implemented for lagrangian reference FEs"
    @unreachable s
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
