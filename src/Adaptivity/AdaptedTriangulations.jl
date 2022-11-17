
"""
  AdaptedTriangulation

  Triangulation produced from an AdaptedDiscreteModel.
  
  Contains: 

  - adapted_model ~> AdaptedDiscreteModel for the triangulation.
  - trian ~> Triangulation for the background model, i.e get_model(adapted_model).
"""
struct AdaptedTriangulation{Dc,Dp,A<:Triangulation{Dc,Dp},B<:AdaptedDiscreteModel} <: Triangulation{Dc,Dp}
  trian::A
  adapted_model::B

  function AdaptedTriangulation(trian::Triangulation{Dc,Dp},model::AdaptedDiscreteModel{Dc2,Dp}) where {Dc,Dc2,Dp}
    @check !isa(trian,AdaptedTriangulation)
    @check Dc <= Dc2
    A = typeof(trian)
    B = typeof(model)
    return new{Dc,Dp,A,B}(trian,model)
  end
end

function get_adapted_model(t::AdaptedTriangulation)
  return t.adapted_model
end

# Relationships
function is_child(t1::AdaptedTriangulation,t2::Triangulation)
  return is_child(get_adapted_model(t1),get_background_model(t2))
end

function is_child(t1::AdaptedTriangulation,t2::AdaptedTriangulation)
  return is_child(get_adapted_model(t1),get_adapted_model(t2))
end

is_child(t1::Triangulation,t2::AdaptedTriangulation) = false

is_related(t1::AdaptedTriangulation,t2::Triangulation) = is_child(t1,t2)
is_related(t1::Triangulation,t2::AdaptedTriangulation) = is_child(t2,t1)
is_related(t1::AdaptedTriangulation,t2::AdaptedTriangulation) = is_child(t1,t2) || is_child(t2,t1)


# Wrap Triangulation API
function Geometry.get_background_model(t::AdaptedTriangulation)
  get_background_model(t.trian) # === get_model(get_adapted_model(t))
end

function Geometry.get_grid(t::AdaptedTriangulation)
  get_grid(t.trian)
end

function Geometry.get_glue(t::AdaptedTriangulation,::Val{d}) where d
  get_glue(t.trian,Val(d))
end

function Base.view(t::AdaptedTriangulation,ids::AbstractArray)
  v = view(t.trian,ids)
  return AdaptedTriangulation(v,t.adapted_model)
end

# Wrap constructors
function Geometry.Triangulation(
  ::Type{ReferenceFE{d}},model::AdaptedDiscreteModel,filter::AbstractArray) where d
  
  trian = Triangulation(ReferenceFE{d},get_model(model),filter)
  return AdaptedTriangulation(trian,model)
end

function Geometry.Triangulation(
  ::Type{ReferenceFE{d}},model::AdaptedDiscreteModel,labels::FaceLabeling;kwargs...) where d
  trian = Triangulation(ReferenceFE{d},get_model(model),labels;kwargs...)
  return AdaptedTriangulation(trian,model)
end

function Geometry.Triangulation(trian::AdaptedTriangulation,args...;kwargs...)
  return AdaptedTriangulation(Triangulation(trian.trian,args...;kwargs...),trian.model)
end

function Geometry.BoundaryTriangulation(model::AdaptedDiscreteModel,args...;kwargs...)
  trian = BoundaryTriangulation(get_model(model),args...;kwargs...)
  return AdaptedTriangulation(trian,model)
end

function Geometry.SkeletonTriangulation(model::AdaptedDiscreteModel,face_to_mask::AbstractVector{Bool})
  trian = SkeletonTriangulation(get_model(model),face_to_mask)
  return AdaptedTriangulation(trian,model)
end

function Geometry.InterfaceTriangulation(model::AdaptedDiscreteModel, cell_to_is_in::Vector{Bool})
  trian = InterfaceTriangulation(get_model(model),cell_to_is_in)
  return AdaptedTriangulation(trian,model)
end

function Geometry.InterfaceTriangulation(model::AdaptedDiscreteModel,cell_to_inout::AbstractVector{<:Integer})
  trian = InterfaceTriangulation(get_model(model),cell_to_inout)
  return AdaptedTriangulation(trian,model)
end

function Geometry.is_change_possible(strian::AdaptedTriangulation,ttrian::AdaptedTriangulation)
  # A) Both Triangulations are exactly the same
  (strian === ttrian) && (return true)

  # B) Same background model -> Default change of Triangulation
  if (get_background_model(strian) === get_background_model(ttrian))
    return is_change_possible(strian.trian,ttrian.trian)
  end

  # C) Different background model, but same type of Triangulation (Skeleton, BodyFitted, View, ...)
  if typeof(strian.trian) == typeof(ttrian.trian) # Is this too restrictive???
    return is_related(strian,ttrian)
  end

  # D) Different background model AND different type of triangulation
  @notimplemented
  return false
end

function Geometry.is_change_possible(strian::AdaptedTriangulation,ttrian::Triangulation)
  # A) Same background model -> Default change of Triangulation
  if (get_background_model(strian) === get_background_model(ttrian))
    return is_change_possible(strian.trian,ttrian)
  end
  
  # B) Different background model, but same type of Triangulation (Skeleton, BodyFitted, View, ...)
  if typeof(strian.trian) == typeof(ttrian)
    return is_related(strian,ttrian)
  end

  # C) Different background model AND different type of triangulation
  @notimplemented
  return false
end

function Geometry.is_change_possible(strian::Triangulation,ttrian::AdaptedTriangulation)
  # A) Same background model -> Default change of Triangulation
  if (get_background_model(strian) === get_background_model(ttrian))
    return is_change_possible(strian,ttrian.trian)
  end
  
  # B) Different background model, but same type of Triangulation (Skeleton, BodyFitted, View, ...)
  if typeof(strian) == typeof(ttrian.trian)
    return is_related(strian,ttrian)
  end

  # C) Different background model AND different type of triangulation
  @notimplemented
  return false
end

function Geometry.best_target(strian::AdaptedTriangulation,ttrian::AdaptedTriangulation)
  @check is_change_possible(strian,ttrian)

  # A) Both Triangulations are exactly the same
  (strian === ttrian) && (return ttrian)

  # B) Same background model -> Default change of Triangulation
  if (get_background_model(strian) === get_background_model(ttrian))
    return best_target(strian.trian,ttrian.trian)
  end

  # C) Different background model, but same type of Triangulation (Skeleton, BodyFitted, View, ...)
  if typeof(strian.trian) == typeof(ttrian.trian)
    is_child(ttrian,strian) ? (return ttrian) : (return strian)
  end

  # D) Different background model AND different type of triangulation
  @notimplemented
  return nothing
end

function Geometry.best_target(strian::AdaptedTriangulation,ttrian::Triangulation)
  @check is_change_possible(strian,ttrian)
  return strian
end

function Geometry.best_target(strian::Triangulation,ttrian::AdaptedTriangulation)
  @check is_change_possible(strian,ttrian)
  return ttrian
end


"""
  Given a AdaptedTriangulation and a CellField defined on the parent(coarse) mesh, 
  returns an equivalent CellField on the fine mesh.
"""
function change_domain_c2f(f_coarse, ftrian::AdaptedTriangulation{Dc,Dp}) where {Dc,Dp}
  @notimplementedif num_dims(get_triangulation(f_coarse)) != Dc

  model  = get_adapted_model(ftrian)
  glue   = get_adaptivity_glue(model)
  if (num_cells(ftrian) != 0)
    # Coarse field but with fine indexing, i.e 
    #   f_f2c[i_fine] = f_coarse[coarse_parent(i_fine)]
    fcell_to_ccell = glue.f2c_faces_map[Dc+1]
    m = Reindex(get_data(f_coarse))
    f_f2c = lazy_map(m,fcell_to_ccell)

    # Fine to coarse coordinate map: x_coarse = Φ^(-1)(x_fine)
    ref_coord_map = get_f2c_reference_coordinate_map(glue)

    # Final map: f_fine(x_fine) = f_f2c ∘ Φ^(-1)(x_fine) = f_coarse(x_coarse)
    f_fine = lazy_map(∘,f_f2c,ref_coord_map)
    return GenericCellField(f_fine,ftrian,ReferenceDomain())
  else
    f_fine = Fill(Gridap.Fields.ConstantField(0.0),num_cells(ftrian))
    return GenericCellField(f_fine,ftrian,ReferenceDomain())
  end
end

function change_domain_f2c(f_fine, ftrian::AdaptedTriangulation{Dc,Dp}) where {Dc,Dp}
  @notimplemented
end

function CellData.change_domain(a::CellField,strian::Triangulation,::ReferenceDomain,ttrian::AdaptedTriangulation,::ReferenceDomain)
  if strian === ttrian
    return a
  end
  @assert is_change_possible(strian,ttrian)

  if (get_background_model(strian) === get_background_model(ttrian))
    return change_domain(a,strian,ReferenceDomain(),ttrian.trian,ReferenceDomain())
  end
  
  return change_domain_c2f(a,ttrian)
end

function CellData.change_domain(a::CellField,strian::AdaptedTriangulation,::ReferenceDomain,ttrian::Triangulation,::ReferenceDomain)
  if strian === ttrian
    return a
  end
  @assert is_change_possible(strian,ttrian)

  if (get_background_model(strian) === get_background_model(ttrian))
    return change_domain(a,strian.trian,ReferenceDomain(),ttrian,ReferenceDomain())
  end
  
  return change_domain_f2c(a,ttrian)
end

function CellData.change_domain(a::CellField,strian::AdaptedTriangulation,::ReferenceDomain,ttrian::AdaptedTriangulation,::ReferenceDomain)
  if is_child(strian,ttrian) # fine to coarse
    return change_domain(a,strian,ReferenceDomain(),ttrian.trian,ReferenceDomain())
  else # coarse to fine
    return change_domain(a,strian.trian,ReferenceDomain(),ttrian,ReferenceDomain())
  end
end

for sdomain in [:ReferenceDomain,:PhysicalDomain]
  for (stype,ttype) in [(:AdaptedTriangulation,:AdaptedTriangulation),(:AdaptedTriangulation,:Triangulation),(:Triangulation,:AdaptedTriangulation)]
    @eval begin
      function CellData.change_domain(a::CellField,strian::$stype,::$sdomain,ttrian::$ttype,::PhysicalDomain)
        a_ref  = change_domain(a,ReferenceDomain())
        atrian = change_domain(a_ref,strian,ReferenceDomain(),ttrian,ReferenceDomain())
        return change_domain(atrian,PhysicalDomain())
      end
    end
  end
end

function Geometry.move_contributions(scell_to_val::AbstractArray, strian::AdaptedTriangulation, ttrian::Triangulation)
  smodel = get_adapted_model(strian)
  @check get_parent(smodel) === get_background_model(ttrian)

  tcell_to_val = move_contributions(scell_to_val,get_adaptivity_glue(smodel))
  return tcell_to_val
end

function Geometry.move_contributions(scell_to_val::AbstractArray, glue::AdaptivityGlue)
  tcell_to_scells = glue.c2f_faces_map
  k = Geometry.CombineContributionsMap(scell_to_val)
  tcell_to_val = lazy_map(k,tcell_to_scells)
  return tcell_to_val
end