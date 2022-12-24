"""
"""
struct CellDof{DS} <: CellDatum
  cell_dof::AbstractArray{<:Union{Dof,AbstractArray{<:Dof}}}
  trian::Triangulation
  domain_style::DS
end

get_data(f::CellDof) = f.cell_dof
get_triangulation(f::CellDof) = f.trian
DomainStyle(::Type{CellDof{DS}}) where DS = DS()
Base.copy(f::CellDof{DS}) where DS = CellDof{DS}(copy(f.cell_dof), f.trian, f.domain_style)

function change_domain(a::CellDof,::ReferenceDomain,::PhysicalDomain)
  @notimplemented
end

function change_domain(a::CellDof,::PhysicalDomain,::ReferenceDomain)
  @notimplemented
end

function change_domain(a::CellDof,ttrian::Triangulation,target_domain::DomainStyle)
  change_domain(a,DomainStyle(a),ttrian,target_domain)
end

function change_domain(a::CellDof,::ReferenceDomain,ttrian::Triangulation,::ReferenceDomain)
  msg = """\n
  We cannot move the given CellDof to the reference domain of the requested triangulation.
  Make sure that the given triangulation is either the same as the triangulation on which the
  CellDof is defined, or that the latter triangulation is the background of the former.
  """
  strian = get_triangulation(a)
  if strian === ttrian
    return a
  end
  @assert is_change_possible(strian,ttrian) msg
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  change_domain_ref_ref(a,ttrian,sglue,tglue)
end

function change_domain(a::CellDof,::PhysicalDomain,ttrian::Triangulation,::PhysicalDomain)
  msg = """\n
  We cannot move the given CellDof to the physical domain of the requested triangulation.
  Make sure that the given triangulation is either the same as the triangulation on which the
  CellDof is defined, or that the latter triangulation is the background of the former.
  """
  strian = get_triangulation(a)
  if strian === ttrian
    return a
  end
  @assert is_change_possible(strian,ttrian) msg
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  change_domain_phys_phys(a,ttrian,sglue,tglue)
end

function change_domain(a::CellDof,::PhysicalDomain,trian::Triangulation,::ReferenceDomain)
  a_trian = change_domain(a,trian,PhysicalDomain())
  change_domain(a_trian,ReferenceDomain())
end

function change_domain(a::CellDof,::ReferenceDomain,trian::Triangulation,::PhysicalDomain)
  a_phys = change_domain(a,PhysicalDomain())
  change_domain(a_phys,trian,PhysicalDomain())
end

function change_domain_ref_ref(
  a::CellDof,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::FaceToFaceGlue)
  sface_to_dof = get_data(a)
  mface_to_sface = sglue.mface_to_tface
  tface_to_mface = tglue.tface_to_mface
  tface_to_mface_map = tglue.tface_to_mface_map
  mface_to_dof = extend(sface_to_dof,mface_to_sface)
  tface_to_dof = lazy_map(Reindex(mface_to_dof),tface_to_mface)
  @notimplementedif ! isa(tface_to_mface_map,AbstractArray{<:GenericField{typeof(identity)}})
  CellDof(tface_to_dof,ttrian,ReferenceDomain())
end

function change_domain_phys_phys(
  a::CellDof,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::FaceToFaceGlue)
  sface_to_dof = get_data(a)
  mface_to_sface = sglue.mface_to_tface
  tface_to_mface = tglue.tface_to_mface
  mface_to_dof = extend(sface_to_dof,mface_to_sface)
  tface_to_dof = lazy_map(Reindex(mface_to_dof),tface_to_mface)
  CellDof(tface_to_dof,ttrian,PhysicalDomain())
end

function get_cell_points(dofs::CellDof)
  cell_to_x = lazy_map(get_nodes,get_data(dofs))
  CellPoint(cell_to_x,get_triangulation(dofs),DomainStyle(dofs))
end

# Evaluation of CellDof

(a::CellDof)(f) = evaluate(a,f)

function evaluate!(cache,s::CellDof,f::CellField)

  trian_f = get_triangulation(f)
  trian_s = get_triangulation(s)

  if trian_f !== trian_s
    @unreachable """\n
    A CellDof can only be evaluated on a CellField defined on the same Triangulation.
    """
  end

  b = change_domain(f,s.domain_style)
  lazy_map(evaluate,get_data(s),get_data(b))
end

function evaluate!(cache, ::CellField, ::CellDof)
  @unreachable """\n
  CellField (f) objects cannot be evaluated at CellDof (s) objects.
  However, CellDofs objects can be evaluated at CellField objects.
  Did you mean evaluate(f,s) instead of evaluate(s,f), i.e.
  f(s) instead of s(f)?
  """
end
