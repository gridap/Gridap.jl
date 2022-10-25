

struct RefinedTriangulation{Dc,Dp,A<:Triangulation{Dc,Dp},B<:RefinedDiscreteModel} <: Triangulation{Dc,Dp}
  trian::A
  refined_model::B

  function RefinedTriangulation(trian::Triangulation{Dc,Dp},model::RefinedDiscreteModel{Dc2,Dp}) where {Dc,Dc2,Dp}
    @check !isa(trian,RefinedTriangulation)
    @check Dc <= Dc2
    A = typeof(trian)
    B = typeof(model)
    return new{Dc,Dp,A,B}(trian,model)
  end
end

function get_refined_model(t::RefinedTriangulation)
  return t.refined_model
end

# Wrap Triangulation API
function Geometry.get_background_model(t::RefinedTriangulation)
  get_background_model(t.trian)
end

function Geometry.get_grid(t::RefinedTriangulation)
  get_grid(t.trian)
end

function Geometry.get_glue(t::RefinedTriangulation,::Val{d}) where d
  get_glue(t.trian,Val(d))
end

function Base.view(t::RefinedTriangulation,ids::AbstractArray)
  v = view(t.trian,ids)
  return RefinedTriangulation(v,t.model)
end

# Wrap constructors
function Geometry.Triangulation(
  ::Type{ReferenceFE{d}},model::RefinedDiscreteModel,filter::AbstractArray) where d
  
  trian = Triangulation(ReferenceFE{d},get_model(model),filter)
  return RefinedTriangulation(trian,model)
end

function Geometry.Triangulation(
  ::Type{ReferenceFE{d}},model::RefinedDiscreteModel,labels::FaceLabeling;kwargs...) where d
  trian = Triangulation(ReferenceFE{d},get_model(model),labels;kwargs...)
  return RefinedTriangulation(trian,model)
end

function Geometry.Triangulation(trian::RefinedTriangulation,args...;kwargs...)
  return RefinedTriangulation(Triangulation(trian.trian,args...;kwargs...),trian.model)
end


function Geometry.is_change_possible(strian::RefinedTriangulation,ttrian::RefinedTriangulation)
  # A) Both Triangulations are exactly the same
  (strian === ttrian) && (return true)

  # B) Same background model -> Default change of Triangulation
  if (get_background_model(strian) === get_background_model(ttrian))
    return is_change_possible(strian.trian,ttrian.trian)
  end

  # C) Different background model, but same type of Triangulation (Skeleton, BodyFitted, View, ...)
  if typeof(strian.trian) == typeof(ttrian.trian) # Is this too restrictive???
    smodel = get_refined_model(strian)
    tmodel = get_refined_model(ttrian)
    a = get_parent(tmodel) === get_model(smodel) # tmodel = refine(smodel)
    b = get_parent(smodel) === get_model(tmodel) # smodel = refine(tmodel)
    return a || b
  end

  # D) Different background model AND different type of triangulation
  @notimplemented
  return false
end

function Geometry.is_change_possible(strian::RefinedTriangulation,ttrian::Triangulation)
  # A) Same background model -> Default change of Triangulation
  if (get_background_model(strian) === get_background_model(ttrian))
    return is_change_possible(strian.trian,ttrian)
  end
  
  # B) Different background model, but same type of Triangulation (Skeleton, BodyFitted, View, ...)
  if typeof(strian.trian) == typeof(ttrian)
    smodel = get_refined_model(strian)
    tmodel = get_background_model(ttrian)
    return get_parent(smodel) === tmodel # smodel = refine(tmodel)
  end

  # C) Different background model AND different type of triangulation
  @notimplemented
  return false
end

# TODO: Are we sure this is symmetric?
function Geometry.is_change_possible(strian::Triangulation,ttrian::RefinedTriangulation)
  return is_change_possible(ttrian,strian)
end

function Geometry.best_target(strian::RefinedTriangulation,ttrian::RefinedTriangulation)
  @check is_change_possible(strian,ttrian)

  # A) Both Triangulations are exactly the same
  (strian === ttrian) && (return ttrian)

  # B) Same background model -> Default change of Triangulation
  if (get_background_model(strian) === get_background_model(ttrian))
    return best_target(strian,ttrian)
  end

  # C) Different background model, but same type of Triangulation (Skeleton, BodyFitted, View, ...)
  if typeof(strian.trian) == typeof(ttrian.trian)
    smodel = get_refined_model(strian)
    tmodel = get_refined_model(ttrian)
    a = get_parent(tmodel) === get_model(smodel) # tmodel = refine(smodel)
    a ? (return ttrian) : (return strian)
  end

  # D) Different background model AND different type of triangulation
  @notimplemented
  return nothing
end

function Geometry.best_target(strian::RefinedTriangulation,ttrian::Triangulation)
  @check is_change_possible(strian,ttrian)
  return strian
end

function Geometry.best_target(strian::Triangulation,ttrian::RefinedTriangulation)
  @check is_change_possible(strian,ttrian)
  return ttrian
end


"""
  Given a RefinedTriangulation and a CellField defined on the parent(coarse) mesh, 
  returns an equivalent CellField on the fine mesh.
"""
function change_domain_c2f(f_coarse, ftrian::RefinedTriangulation{Dc,Dp}) where {Dc,Dp}
  @notimplementedif num_dims(get_triangulation(f_coarse)) != Dc

  model  = get_refined_model(ftrian)
  glue   = get_refinement_glue(model)
  if (num_cells(ftrian) != 0)
    # Coarse field but with fine indexing, i.e 
    #   f_f2c[i_fine] = f_coarse[coarse_parent(i_fine)]
    fcell_to_ccell = glue.f2c_faces_map[Dc+1]
    m = Reindex(get_data(f_coarse))
    f_f2c = lazy_map(m,fcell_to_ccell)

    # Fine to coarse coordinate map: x_coarse = Φ^(-1)(x_fine)
    ref_coord_map = get_f2c_ref_coordinate_map(glue)

    # Final map: f_fine(x_fine) = f_f2c ∘ Φ^(-1)(x_fine) = f_coarse(x_coarse)
    f_fine = lazy_map(∘,f_f2c,ref_coord_map)
    return GenericCellField(f_fine,ftrian,ReferenceDomain())
  else
    f_fine = Fill(Gridap.Fields.ConstantField(0.0),num_cells(ftrian))
    return GenericCellField(f_fine,ftrian,ReferenceDomain())
  end
end

function CellData.change_domain(a::CellField,ttrian::RefinedTriangulation,::ReferenceDomain)
  strian = get_triangulation(a)
  if strian === ttrian
    return a
  end
  @assert is_change_possible(strian,ttrian)

  if (get_background_model(strian) === get_background_model(ttrian))
    return change_domain(a,ttrian.trian,ReferenceDomain())
  end
  
  return change_domain_c2f(a,ttrian)
end

function CellData.change_domain(a::CellData.OperationCellField,ttrian::RefinedTriangulation,::ReferenceDomain)
  strian = get_triangulation(a)
  if strian === ttrian
    return a
  end
  @assert is_change_possible(strian,ttrian)

  if (get_background_model(strian) === get_background_model(ttrian))
    return change_domain(a,ttrian.trian,ReferenceDomain())
  end
  
  return change_domain_c2f(a,ttrian)
end

function CellData.change_domain(a::CellField,ttrian::RefinedTriangulation,::PhysicalDomain)
  strian = get_triangulation(a)
  if strian === ttrian
    return a
  end
  @assert is_change_possible(strian,ttrian)

  if (get_background_model(strian) === get_background_model(ttrian))
    return change_domain(a,ttrian.trian,PhysicalDomain())
  end

  @notimplemented
end
