
"""

  Triangulation produced from an AdaptedDiscreteModel.
  
  Contains: 

  - adapted_model :: `AdaptedDiscreteModel` for the triangulation.
  - trian :: `Triangulation` extracted from the background model, i.e `get_model(adapted_model)`.
"""
struct AdaptedTriangulation{Dc,Dp,A<:Triangulation{Dc,Dp},B<:AdaptedDiscreteModel} <: Triangulation{Dc,Dp}
  trian::A
  adapted_model::B

  function AdaptedTriangulation(trian::Triangulation{Dc,Dp},model::AdaptedDiscreteModel{Dc2,Dp}) where {Dc,Dc2,Dp}
    @check !isa(trian,AdaptedTriangulation)
    @check get_background_model(trian) === get_model(model)
    @check Dc <= Dc2
    A = typeof(trian)
    B = typeof(model)
    return new{Dc,Dp,A,B}(trian,model)
  end
end

function get_adapted_model(t::AdaptedTriangulation)
  return t.adapted_model
end

function get_adaptivity_glue(t::AdaptedTriangulation)
  return get_adaptivity_glue(get_adapted_model(t))
end

# Relationships
function is_child(t1::AdaptedTriangulation,t2::Triangulation)
  return is_child(get_adapted_model(t1),get_background_model(t2))
end

function is_child(t1::AdaptedTriangulation,t2::AdaptedTriangulation)
  return is_child(get_adapted_model(t1),get_adapted_model(t2))
end

is_child(t1::Triangulation,t2::AdaptedTriangulation) = false

is_related(t1::Triangulation,t2::Triangulation) = is_child(t1,t2) || is_child(t2,t1)


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

function Geometry.get_facet_normal(trian::AdaptedTriangulation)
  get_facet_normal(trian.trian)
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
  return AdaptedTriangulation(Triangulation(trian.trian,args...;kwargs...),trian.adapted_model)
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

for (stype,ttype) in [(:AdaptedTriangulation,:AdaptedTriangulation),(:AdaptedTriangulation,:Triangulation),(:Triangulation,:AdaptedTriangulation)]
  sstrian = (stype==:AdaptedTriangulation) ? :(strian.trian) : :(strian)
  tttrian = (ttype==:AdaptedTriangulation) ? :(ttrian.trian) : :(ttrian)
  @eval begin
    function CellData.is_change_possible(strian::$stype,ttrian::$ttype)
      (strian === ttrian) && (return true)

      if (get_background_model(strian) === get_background_model(ttrian))
        return is_change_possible($sstrian,$tttrian)
      end
      
      strian_is_cell_wise = (num_cell_dims(strian) == num_point_dims(strian))
      trians_are_related  = is_related(strian,ttrian)
      return strian_is_cell_wise && trians_are_related
    end
  end
  @eval begin
    function Geometry.best_target(strian::$stype,ttrian::$ttype)
      @check is_change_possible(strian,ttrian)
    
      (strian === ttrian) && (return ttrian)
    
      if (get_background_model(strian) === get_background_model(ttrian))
        return best_target($sstrian,$tttrian)
      end

      ttrian_is_cell_wise = (num_cell_dims(ttrian) == num_point_dims(ttrian))
      strian_is_cell_wise = (num_cell_dims(strian) == num_point_dims(strian))
    
      if strian_is_cell_wise && ttrian_is_cell_wise
        # Choose the child mesh, both ways are possible
        is_child(ttrian,strian) ? (return ttrian) : (return strian)
      elseif strian_is_cell_wise
        return strian
      end

      @notimplemented
      return nothing
    end
  end
end

function CellData.change_domain(a::CellField,strian::Triangulation,::ReferenceDomain,ttrian::AdaptedTriangulation,::ReferenceDomain)
  @check is_change_possible(strian,ttrian)

  if (get_background_model(strian) === get_background_model(ttrian))
    b = change_domain(a,strian,ReferenceDomain(),ttrian.trian,ReferenceDomain())
    return CellData.similar_cell_field(b,get_data(b),ttrian,ReferenceDomain())
  end
  
  return change_domain_o2n(a,strian,ttrian)
end

function CellData.change_domain(a::CellField,strian::AdaptedTriangulation,::ReferenceDomain,ttrian::Triangulation,::ReferenceDomain)
  @check is_change_possible(strian,ttrian)

  if (get_background_model(strian) === get_background_model(ttrian))
    return change_domain(a,strian.trian,ReferenceDomain(),ttrian,ReferenceDomain())
  end
  
  return change_domain_n2o(a,strian,ttrian)
end

function CellData.change_domain(a::CellField,strian::AdaptedTriangulation,::ReferenceDomain,ttrian::AdaptedTriangulation,::ReferenceDomain)
  if strian === ttrian
    return a
  end

  if is_child(strian,ttrian) # fine to coarse
    b = change_domain(a,strian,ReferenceDomain(),ttrian.trian,ReferenceDomain())
    return CellData.similar_cell_field(b,get_data(b),ttrian,ReferenceDomain())
  else # coarse to fine
    return change_domain(a,strian.trian,ReferenceDomain(),ttrian,ReferenceDomain())
  end
end

for sdomain in [:ReferenceDomain,:PhysicalDomain]
  for (stype,ttype) in [(:AdaptedTriangulation,:AdaptedTriangulation),(:AdaptedTriangulation,:Triangulation),(:Triangulation,:AdaptedTriangulation)]
    sstrian = (stype==:AdaptedTriangulation) ? :(strian.trian) : :(strian)
    tttrian = (ttype==:AdaptedTriangulation) ? :(ttrian.trian) : :(ttrian)
    @eval begin
      function CellData.change_domain(a::CellField,strian::$stype,sd::$sdomain,ttrian::$ttype,::PhysicalDomain)
        if (get_background_model(strian) === get_background_model(ttrian))
          return change_domain(a,$sstrian,sd,$tttrian,PhysicalDomain())
        end

        a_ref  = change_domain(a,ReferenceDomain())
        atrian = change_domain(a_ref,strian,ReferenceDomain(),ttrian,ReferenceDomain())
        return change_domain(atrian,PhysicalDomain())
      end
    end
  end
end

function Geometry.move_contributions(scell_to_val::AbstractArray, strian::AdaptedTriangulation, ttrian::Triangulation{Dc}) where Dc
  # Default Gridap 
  if get_background_model(strian) === get_background_model(ttrian)
    return move_contributions(scell_to_val,strian.trian,ttrian)
  end
  
  # Else 
  smodel = get_adapted_model(strian)
  @check get_parent(smodel) === get_background_model(ttrian)

  return move_contributions_f2c(scell_to_val,strian,ttrian)
end

function move_contributions_f2c(fine_tface_to_val::AbstractArray, ftrian::AdaptedTriangulation, ctrian::Triangulation{Dc}) where Dc
  glue  = get_adaptivity_glue(ftrian)
  fglue = get_glue(ftrian,Val(Dc))
  cglue = get_glue(ctrian,Val(Dc))

  # Fine Triangulation -> Fine Model
  fmodel = get_background_model(ftrian)
  cell_ftrian = Triangulation(ReferenceFE{Dc},fmodel)
  fine_mface_to_val = move_contributions(fine_tface_to_val,ftrian.trian,cell_ftrian)

  # Fine model -> Coarse Model
  coarse_mface_to_val = move_contributions(fine_mface_to_val,glue)

  # Coarse Model -> Coarse Triangulation
  coarse_tface_to_val = lazy_map(Reindex(coarse_mface_to_val),cglue.tface_to_mface)

  return coarse_tface_to_val
end

function Geometry.move_contributions(scell_to_val::AbstractArray, glue::AdaptivityGlue)
  tcell_to_scells = glue.o2n_faces_map
  k = Geometry.CombineContributionsMap(scell_to_val)
  tcell_to_val = lazy_map(k,tcell_to_scells)
  return tcell_to_val
end

# Change domain o2n/n2o

function change_domain_o2n(f_old,old_trian::Triangulation,new_trian::AdaptedTriangulation)
  glue = get_adaptivity_glue(new_trian)
  return change_domain_o2n(f_old,old_trian,new_trian,glue)
end

function change_domain_n2o(f_new,new_trian::AdaptedTriangulation,old_trian::Triangulation)
  glue = get_adaptivity_glue(new_trian)
  return change_domain_n2o(f_new,new_trian,old_trian,glue)
end

"""
  Given a AdaptivityGlue and a CellField defined on the parent(old) mesh, 
  returns an equivalent CellField on the child(new) mesh.
"""
function change_domain_o2n(f_old,old_trian::Triangulation,new_trian::AdaptedTriangulation,glue::AdaptivityGlue)
  @notimplemented
end

function change_domain_o2n(f_coarse,ctrian::Triangulation{Dc},ftrian::AdaptedTriangulation,glue::AdaptivityGlue{<:RefinementGlue}) where Dc
  cglue = get_glue(ctrian,Val(Dc))
  fglue = get_glue(ftrian,Val(Dc))

  @notimplementedif num_point_dims(ctrian) != Dc
  @notimplementedif isa(fglue,Nothing)

  if (num_cells(ctrian) != 0)
    ### Old Triangulation -> Old Model
    coarse_tface_to_field = CellData.get_data(f_coarse)
    coarse_mface_to_field = extend(coarse_tface_to_field,cglue.mface_to_tface)

    ### Old Model -> New Model
    # Coarse field but with fine indexing, i.e 
    #   f_f2c[i_fine] = f_coarse[coarse_parent(i_fine)]
    f_f2c = o2n_reindex(coarse_mface_to_field,glue)

    # Fine to coarse coordinate map: x_coarse = Φ^(-1)(x_fine)
    ref_coord_map = get_n2o_reference_coordinate_map(glue)

    # f_fine(x_fine) = f_f2c ∘ Φ^(-1)(x_fine) = f_coarse(x_coarse)
    fine_mface_to_field = lazy_map(Broadcasting(∘),f_f2c,ref_coord_map)

    ### New Model -> New Triangulation
    fine_tface_to_field = lazy_map(Reindex(fine_mface_to_field),fglue.tface_to_mface)
    f_fine = lazy_map(Broadcasting(∘),fine_tface_to_field,fglue.tface_to_mface_map)

    return CellData.similar_cell_field(f_coarse,f_fine,ftrian,ReferenceDomain())
  else
    f_fine = Fill(Fields.ConstantField(0.0),num_cells(ftrian))
    return CellData.similar_cell_field(f_coarse,f_fine,ftrian,ReferenceDomain())
  end
end

function change_domain_o2n(
  f_old,
  old_trian::Triangulation{Dc},
  new_trian::AdaptedTriangulation,
  glue::AdaptivityGlue{<:MixedGlue}) where {Dc}

  oglue = get_glue(old_trian,Val(Dc))
  nglue = get_glue(new_trian,Val(Dc))

  @notimplementedif num_point_dims(old_trian) != Dc
  @notimplementedif isa(nglue,Nothing)

  if (num_cells(old_trian) != 0)
    # If mixed refinement/coarsening, then f_c2f is a Table
    f_old_data  = CellData.get_data(f_old)
    f_c2f       = o2n_reindex(f_old_data,glue)
    new_rrules  = get_new_cell_refinement_rules(glue)
    field_array = lazy_map(OldToNewField, f_c2f, new_rrules, glue.n2o_cell_to_child_id)
    return CellData.similar_cell_field(f_old,field_array,new_trian,ReferenceDomain())
  else
    f_new = Fill(Fields.ConstantField(0.0),num_cells(new_trian))
    return CellData.similar_cell_field(f_old,f_new,new_trian,ReferenceDomain())
  end 
end

"""
  Given a AdaptivityGlue and a CellField defined on the child(new) mesh, 
  returns an equivalent CellField on the parent(old) mesh.
"""
function change_domain_n2o(f_new,new_trian::AdaptedTriangulation,old_trian::Triangulation,glue::AdaptivityGlue)
  @notimplemented
end

function change_domain_n2o(f_fine,ftrian::AdaptedTriangulation{Dc},ctrian::Triangulation,glue::AdaptivityGlue{<:RefinementGlue}) where Dc
  fglue = get_glue(ftrian,Val(Dc))
  cglue = get_glue(ctrian,Val(Dc))

  @notimplementedif num_point_dims(ftrian) != Dc
  @notimplementedif isa(cglue,Nothing)

  if (num_cells(ctrian) != 0)
    ### New Triangulation -> New Model
    fine_tface_to_field = CellData.get_data(f_fine)
    fine_mface_to_field = extend(fine_tface_to_field,fglue.mface_to_tface)

    ### New Model -> Old Model
    # f_c2f[i_coarse] = [f_fine[i_fine_1], ..., f_fine[i_fine_nChildren]]
    f_c2f = n2o_reindex(fine_mface_to_field,glue)

    child_ids = n2o_reindex(glue.n2o_cell_to_child_id,glue)
    rrules    = get_old_cell_refinement_rules(glue)
    coarse_mface_to_field = lazy_map(FineToCoarseField,f_c2f,rrules,child_ids)

    ### Old Model -> Old Triangulation
    coarse_tface_to_field = lazy_map(Reindex(coarse_mface_to_field),cglue.tface_to_mface)
    f_coarse = lazy_map(Broadcasting(∘),coarse_tface_to_field,cglue.tface_to_mface_map)

    return CellData.similar_cell_field(f_fine,f_coarse,ctrian,ReferenceDomain())
  else
    f_coarse = Fill(Fields.ConstantField(0.0),num_cells(fcoarse))
    return CellData.similar_cell_field(f_fine,f_coarse,ctrian,ReferenceDomain())
  end
end


# Specialisation for Skeleton Pairs

function change_domain_o2n(f_old,old_trian::Triangulation,new_trian::AdaptedTriangulation{Dc,Dp,<:SkeletonTriangulation}) where {Dc,Dp}
  @check isa(CellData.get_data(f_old),Fill)
  new_trian_plus = AdaptedTriangulation(new_trian.trian.plus,new_trian.adapted_model)
  return change_domain_o2n(f_old,old_trian,new_trian_plus)
end

function change_domain_n2o(f_new,new_trian::AdaptedTriangulation,old_trian::SkeletonTriangulation)
  @check isa(CellData.get_data(f_new),Fill)
  return change_domain_n2o(f_new,new_trian,old_trian.plus)
end

function change_domain_o2n(f_old::CellData.CellFieldAt,old_trian::Triangulation,new_trian::AdaptedTriangulation{Dc,Dp,<:SkeletonTriangulation}) where {Dc,Dp}
  if isa(f_old,CellData.CellFieldAt{:plus})
    _new_trian = AdaptedTriangulation(new_trian.trian.plus,new_trian.adapted_model)
  elseif isa(f_old,CellData.CellFieldAt{:minus})
    _new_trian = AdaptedTriangulation(new_trian.trian.minus,new_trian.adapted_model)
  else
    @unreachable
  end
  return change_domain_o2n(f_old,old_trian,_new_trian)
end

function change_domain_n2o(f_new::CellData.CellFieldAt,new_trian::AdaptedTriangulation,old_trian::SkeletonTriangulation)
  if isa(f_new,CellData.CellFieldAt{:plus})
    _old_trian = old_trian.plus
  elseif isa(f_new,CellData.CellFieldAt{:minus})
    _old_trian = old_trian.minus
  else
    @unreachable
  end
  return change_domain_n2o(f_new,new_trian,_old_trian)
end


