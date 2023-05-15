

struct GluedTriangulation{Dc,Dp,A,B,C} <: Triangulation{Dc,Dp}
  trian::A
  parent_trian::B
  glue::C

  function GluedTriangulation(trian::Triangulation{Dc,Dp},
                              parent_trian::Triangulation{Dc},
                              glue::AdaptivityGlue) where {Dc,Dp}
    
    A = typeof(trian)
    B = typeof(parent_trian)
    C = typeof(glue)
    return new{Dc,Dp,A,B,C}(trian,parent_trian,glue)
  end
end

function get_adaptivity_glue(t::GluedTriangulation)
  return t.glue
end

# Relationships
function is_child(t1::GluedTriangulation,t2::Triangulation)
  return t2 === t1.parent_trian
end

function is_child(t1::GluedTriangulation,t2::GluedTriangulation)
  return is_child(t1,t2.trian)
end

is_child(t1::Triangulation,t2::GluedTriangulation) = false

# Wrap Triangulation API
function Geometry.get_background_model(t::GluedTriangulation)
  get_background_model(t.trian) # === get_model(get_adapted_model(t))
end

function Geometry.get_grid(t::GluedTriangulation)
  get_grid(t.trian)
end

function Geometry.get_glue(t::GluedTriangulation,::Val{d}) where d
  get_glue(t.trian,Val(d))
end

function Geometry.get_facet_normal(trian::GluedTriangulation)
  get_facet_normal(trian.trian)
end

for (stype,ttype) in [(:GluedTriangulation,:GluedTriangulation),(:GluedTriangulation,:Triangulation),(:Triangulation,:GluedTriangulation)]
  sstrian = (stype==:GluedTriangulation) ? :(strian.trian) : :(strian)
  tttrian = (ttype==:GluedTriangulation) ? :(ttrian.trian) : :(ttrian)
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

function CellData.change_domain(a::CellField,strian::Triangulation,::ReferenceDomain,ttrian::GluedTriangulation,::ReferenceDomain)
  @check is_change_possible(strian,ttrian)

  if (get_background_model(strian) === get_background_model(ttrian))
    b = change_domain(a,strian,ReferenceDomain(),ttrian.trian,ReferenceDomain())
    return CellData.similar_cell_field(b,get_data(b),ttrian,ReferenceDomain())
  end
  
  return change_domain_o2n(a,strian,ttrian)
end

function CellData.change_domain(a::CellField,strian::GluedTriangulation,::ReferenceDomain,ttrian::Triangulation,::ReferenceDomain)
  @check is_change_possible(strian,ttrian)

  if (get_background_model(strian) === get_background_model(ttrian))
    return change_domain(a,strian.trian,ReferenceDomain(),ttrian,ReferenceDomain())
  end
  
  return change_domain_n2o(a,strian,ttrian)
end

function CellData.change_domain(a::CellField,strian::GluedTriangulation,::ReferenceDomain,ttrian::GluedTriangulation,::ReferenceDomain)
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
  for (stype,ttype) in [(:GluedTriangulation,:GluedTriangulation),(:GluedTriangulation,:Triangulation),(:Triangulation,:GluedTriangulation)]
    @eval begin
      function CellData.change_domain(a::CellField,strian::$stype,::$sdomain,ttrian::$ttype,::PhysicalDomain)
        a_ref  = change_domain(a,ReferenceDomain())
        atrian = change_domain(a_ref,strian,ReferenceDomain(),ttrian,ReferenceDomain())
        return change_domain(atrian,PhysicalDomain())
      end
    end
  end
end

function Geometry.move_contributions(scell_to_val::AbstractArray, strian::GluedTriangulation, ttrian::Triangulation{Dc}) where Dc
  # Default Gridap 
  if get_background_model(strian) === get_background_model(ttrian)
    return move_contributions(scell_to_val,strian.trian,ttrian)
  end
  
  # Else 
  smodel = get_adapted_model(strian)
  @check get_parent(smodel) === get_background_model(ttrian)

  return move_contributions_f2c(scell_to_val,strian,ttrian)
end

function move_contributions_f2c(fine_tface_to_val::AbstractArray, ftrian::GluedTriangulation, ctrian::Triangulation{Dc}) where Dc
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

# get_cell_dof_ids 

function FESpaces.get_cell_dof_ids(f::FESpaces.SingleFieldFESpace,ttrian::Adaptivity.GluedTriangulation)
  strian = get_triangulation(f)
  if strian == ttrian.parent_trian
    coarse_dof_ids = get_cell_dof_ids(f)
    return Adaptivity.c2f_reindex(coarse_dof_ids,ttrian.glue)
  else
    return get_cell_dof_ids(f,ttrian.trian)
  end
end

# Change domain o2n/n2o

function change_domain_o2n(f_old,old_trian::Triangulation,new_trian::GluedTriangulation)
  glue = get_adaptivity_glue(new_trian)
  return change_domain_o2n(f_old,old_trian,new_trian,glue)
end

function change_domain_n2o(f_new,new_trian::GluedTriangulation,old_trian::Triangulation)
  glue = get_adaptivity_glue(new_trian)
  return change_domain_n2o(f_new,new_trian,old_trian,glue)
end

"""
  Given a AdaptivityGlue and a CellField defined on the parent(old) mesh, 
  returns an equivalent CellField on the child(new) mesh.
"""
function change_domain_o2n(f_old,old_trian::Triangulation,new_trian::GluedTriangulation,glue::AdaptivityGlue)
  @notimplemented
end

function change_domain_o2n(f_coarse,ctrian::Triangulation{Dc},ftrian::GluedTriangulation,glue::AdaptivityGlue{<:RefinementGlue}) where Dc
  cglue = get_glue(ctrian,Val(Dc))
  fglue = get_glue(ftrian,Val(Dc))

  @notimplementedif num_point_dims(ctrian) != Dc
  @notimplementedif isa(fglue,Nothing)

  if (num_cells(ctrian) != 0)
    ### Old Triangulation -> Old Model
    coarse_tface_to_field = CellData.get_data(f_coarse)
    #coarse_mface_to_field = extend(coarse_tface_to_field,cglue.mface_to_tface)

    ### Old Model -> New Model
    # Coarse field but with fine indexing, i.e 
    #   f_f2c[i_fine] = f_coarse[coarse_parent(i_fine)]
    f_f2c = c2f_reindex(coarse_tface_to_field,glue)

    # Fine to coarse coordinate map: x_coarse = Φ^(-1)(x_fine)
    ref_coord_map = get_n2o_reference_coordinate_map(glue)

    # f_fine(x_fine) = f_f2c ∘ Φ^(-1)(x_fine) = f_coarse(x_coarse)
    f_fine = lazy_map(∘,f_f2c,ref_coord_map)

    ### New Model -> New Triangulation
    #fine_tface_to_field = lazy_map(Reindex(fine_mface_to_field),fglue.tface_to_mface)
    #f_fine = lazy_map(Broadcasting(∘),fine_tface_to_field,fglue.tface_to_mface_map)

    return CellData.similar_cell_field(f_coarse,f_fine,ftrian,ReferenceDomain())
  else
    f_fine = Fill(Gridap.Fields.ConstantField(0.0),num_cells(ftrian))
    return CellData.similar_cell_field(f_coarse,f_fine,ftrian,ReferenceDomain())
  end
end

"""
  Given a AdaptivityGlue and a CellField defined on the child(new) mesh, 
  returns an equivalent CellField on the parent(old) mesh.
"""
function change_domain_n2o(f_new,new_trian::GluedTriangulation,old_trian::Triangulation,glue::AdaptivityGlue)
  @notimplemented
end

function change_domain_n2o(f_fine,ftrian::GluedTriangulation{Dc},ctrian::Triangulation,glue::AdaptivityGlue{<:RefinementGlue}) where Dc
  fglue = get_glue(ftrian,Val(Dc))
  cglue = get_glue(ctrian,Val(Dc))

  @notimplementedif num_point_dims(ftrian) != Dc
  @notimplementedif isa(cglue,Nothing)

  if (num_cells(ctrian) != 0)
    ### New Triangulation -> New Model
    fine_tface_to_field = CellData.get_data(f_fine)
    #fine_mface_to_field = extend(fine_tface_to_field,fglue.mface_to_tface)

    ### New Model -> Old Model
    # f_c2f[i_coarse] = [f_fine[i_fine_1], ..., f_fine[i_fine_nChildren]]
    f_c2f = f2c_reindex(fine_tface_to_field,glue)

    rrules   = get_old_cell_refinement_rules(glue)
    f_coarse = lazy_map(FineToCoarseField,f_c2f,rrules)

    ### Old Model -> Old Triangulation
    #coarse_tface_to_field = lazy_map(Reindex(coarse_mface_to_field),cglue.tface_to_mface)
    #f_coarse = lazy_map(Broadcasting(∘),coarse_tface_to_field,cglue.tface_to_mface_map)

    return CellData.similar_cell_field(f_fine,f_coarse,ctrian,ReferenceDomain())
  else
    f_coarse = Fill(Gridap.Fields.ConstantField(0.0),num_cells(fcoarse))
    return CellData.similar_cell_field(f_fine,f_coarse,ctrian,ReferenceDomain())
  end
end


# Specialisation for Skeleton Pairs
"""
function change_domain_o2n(f_old,old_trian::Triangulation,new_trian::GluedTriangulation{Dc,Dp,<:SkeletonTriangulation}) where {Dc,Dp}
  @check isa(CellData.get_data(f_old),Fill)
  new_trian_plus = GluedTriangulation(new_trian.trian.plus,new_trian.adapted_model)
  return change_domain_o2n(f_old,old_trian,new_trian_plus)
end

function change_domain_n2o(f_new,new_trian::GluedTriangulation,old_trian::SkeletonTriangulation)
  @check isa(CellData.get_data(f_new),Fill)
  return change_domain_n2o(f_new,new_trian,old_trian.plus)
end

function change_domain_o2n(f_old::CellData.CellFieldAt,old_trian::Triangulation,new_trian::GluedTriangulation{Dc,Dp,<:SkeletonTriangulation}) where {Dc,Dp}
  if isa(f_old,CellData.CellFieldAt{:plus})
    _new_trian = GluedTriangulation(new_trian.trian.plus,new_trian.adapted_model)
  elseif isa(f_old,CellData.CellFieldAt{:minus})
    _new_trian = GluedTriangulation(new_trian.trian.minus,new_trian.adapted_model)
  else
    @unreachable
  end
  return change_domain_o2n(f_old,old_trian,_new_trian)
end

function change_domain_n2o(f_new::CellData.CellFieldAt,new_trian::GluedTriangulation,old_trian::SkeletonTriangulation)
  if isa(f_new,CellData.CellFieldAt{:plus})
    _old_trian = old_trian.plus
  elseif isa(f_new,CellData.CellFieldAt{:minus})
    _old_trian = old_trian.minus
  else
    @unreachable
  end
  return change_domain_n2o(f_new,new_trian,_old_trian)
end
"""

