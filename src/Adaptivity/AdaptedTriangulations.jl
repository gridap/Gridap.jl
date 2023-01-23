
"""
  AdaptedTriangulation

  Triangulation produced from an AdaptedDiscreteModel.
  
  Contains: 

  - adapted_model ~> AdaptedDiscreteModel for the triangulation.
  - trian ~> Triangulation extracted from the background model, i.e get_model(adapted_model).
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
      
      if isa($sstrian,BodyFittedTriangulation) && isa($tttrian,BodyFittedTriangulation)
        return is_related(strian,ttrian)
      end

      # Support for Boundary, Skeletton, ... triangulations is not yet implemented
      # unless they are based on the same background model.
      @notimplemented
      return false
    end
  end
  @eval begin
    function Geometry.best_target(strian::$stype,ttrian::$ttype)
      @check is_change_possible(strian,ttrian)
    
      (strian === ttrian) && (return ttrian)
    
      if (get_background_model(strian) === get_background_model(ttrian))
        return best_target($sstrian,$tttrian)
      end
    
      if isa($sstrian,BodyFittedTriangulation) && isa($tttrian,BodyFittedTriangulation)
        is_child(ttrian,strian) ? (return ttrian) : (return strian)
      end
    
      # Support for Boundary, Skeletton, ... triangulations is not yet implemented
      # unless they are based on the same background model.
      @notimplemented
      return nothing
    end
  end
end

function CellData.change_domain(a::CellField,strian::Triangulation,::ReferenceDomain,ttrian::AdaptedTriangulation,::ReferenceDomain)
  if strian === ttrian
    return a
  end
  @check is_change_possible(strian,ttrian)

  if (get_background_model(strian) === get_background_model(ttrian))
    return change_domain(a,strian,ReferenceDomain(),ttrian.trian,ReferenceDomain())
  end
  
  return change_domain_o2n(a,ttrian)
end

function CellData.change_domain(a::CellField,strian::AdaptedTriangulation,::ReferenceDomain,ttrian::Triangulation,::ReferenceDomain)
  if strian === ttrian
    return a
  end
  @check is_change_possible(strian,ttrian)

  if (get_background_model(strian) === get_background_model(ttrian))
    return change_domain(a,strian.trian,ReferenceDomain(),ttrian,ReferenceDomain())
  end
  
  return change_domain_n2o(a,ttrian)
end

function CellData.change_domain(a::CellField,strian::AdaptedTriangulation,::ReferenceDomain,ttrian::AdaptedTriangulation,::ReferenceDomain)
  if is_child(strian,ttrian) # fine to coarse
    b = change_domain(a,strian,ReferenceDomain(),ttrian.trian,ReferenceDomain())
    return CellData.similar_cell_field(b,get_data(b),ttrian,ReferenceDomain())
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
  tcell_to_scells = glue.o2n_faces_map
  k = Geometry.CombineContributionsMap(scell_to_val)
  tcell_to_val = lazy_map(k,tcell_to_scells)
  return tcell_to_val
end

function change_domain_o2n(f_old, new_trian::AdaptedTriangulation{Dc,Dp}) where {Dc,Dp}
  @notimplementedif num_dims(get_triangulation(f_old)) != Dc
  glue = get_adaptivity_glue(new_trian)
  return change_domain_o2n(f_old,new_trian,glue)
end

function change_domain_n2o(f_new, old_trian::Triangulation{Dc,Dp}) where {Dc,Dp}
  new_trian = get_triangulation(f_new)
  @notimplementedif num_dims(new_trian) != Dc
  @check isa(new_trian,AdaptedTriangulation)
  glue = get_adaptivity_glue(new_trian)
  return change_domain_n2o(f_new,old_trian,glue)
end

"""
  Given a AdaptivityGlue and a CellField defined on the parent(old) mesh, 
  returns an equivalent CellField on the child(new) mesh.
"""
function change_domain_o2n(f_old,new_trian::AdaptedTriangulation,glue::AdaptivityGlue)
  @notimplemented
end

function change_domain_o2n(f_coarse,ftrian::AdaptedTriangulation{Dc},glue::AdaptivityGlue{<:RefinementGlue,Dc}) where Dc
  ctrian = get_triangulation(f_coarse)
  @notimplementedif num_dims(ctrian) != Dc

  if (num_cells(ctrian) != 0)
    # Coarse field but with fine indexing, i.e 
    #   f_f2c[i_fine] = f_coarse[coarse_parent(i_fine)]
    f_f2c = c2f_reindex(f_coarse,glue)

    # Fine to coarse coordinate map: x_coarse = Φ^(-1)(x_fine)
    ref_coord_map = get_n2o_reference_coordinate_map(glue)

    # Final map: f_fine(x_fine) = f_f2c ∘ Φ^(-1)(x_fine) = f_coarse(x_coarse)
    f_fine = lazy_map(∘,f_f2c,ref_coord_map)
    return GenericCellField(f_fine,ftrian,ReferenceDomain())
  else
    f_fine = Fill(Gridap.Fields.ConstantField(0.0),num_cells(ftrian))
    return GenericCellField(f_fine,ftrian,ReferenceDomain())
  end
end

"""
  Given a AdaptivityGlue and a CellField defined on the child(new) mesh, 
  returns an equivalent CellField on the parent(old) mesh.
"""
function change_domain_n2o(f_new,old_trian::Triangulation,glue::AdaptivityGlue)
  @notimplemented
end

function change_domain_n2o(f_fine,ctrian::Triangulation{Dc},glue::AdaptivityGlue{<:RefinementGlue,Dc}) where Dc
  @notimplementedif num_dims(ctrian) != Dc
  msg = "Evaluating a fine CellField in the coarse mesh is costly! If you are using this feature 
         to integrate, consider using a CompositeMeasure instead (see test/AdaptivityTests/GridTransferTests.jl)."
  @warn msg

  if (num_cells(ctrian) != 0)
    # f_c2f[i_coarse] = [f_fine[i_fine_1], ..., f_fine[i_fine_nChildren]]
    f_c2f = f2c_reindex(f_fine,glue)

    rrules   = get_old_cell_refinement_rules(glue)
    f_coarse = lazy_map(FineToCoarseField,f_c2f,rrules)
    return GenericCellField(f_coarse,ctrian,ReferenceDomain())
  else
    f_coarse = Fill(Gridap.Fields.ConstantField(0.0),num_cells(ftrian))
    return GenericCellField(f_coarse,ctrian,ReferenceDomain())
  end
end
