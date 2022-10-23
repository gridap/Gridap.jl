

struct RefinedTriangulation{Dc,Dp,A<:Triangulation{Dc,Dp},B<:RefinedDiscreteModel} <: Triangulation{Dc,Dp}
  trian::A
  model::B

  function RefinedTriangulation(trian::Triangulation{Dc,Dp},model::RefinedDiscreteModel{Dc2,Dp}) where {Dc,Dc2,Dp}
    @check Dc <= Dc2
    A = typeof(trian)
    B = typeof(model)
    return new{Dc,Dp,A,B}(trian,model)
  end
end

function get_refined_model(t::RefinedTriangulation)
  return t.model
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

# Domain changes
# TODO: This assumes we have the same type of triangulation on both refinement levels!
#       we might want to change this in the future when doing hybrid methods etc...

function Geometry.is_change_possible(strian::RefinedTriangulation,ttrian::RefinedTriangulation)
  (strian === ttrian) && (return true)
  smodel = get_refined_model(strian)
  tmodel = get_refined_model(ttrian)
  a = get_parent(tmodel) === smodel # tmodel = refine(smodel)
  b = get_parent(smodel) === tmodel # smodel = refine(tmodel)
  return a || b
end

function Geometry.is_change_possible(strian::RefinedTriangulation,ttrian::T) where {T <: Triangulation}
  smodel = get_refined_model(strian)
  tmodel = get_background_model(ttrian)
  a = get_model(smodel) === tmodel  # It is fundamentally the same model
  b = get_parent(smodel) === tmodel # smodel = refine(tmodel)
  return a || b
end

function Geometry.is_change_possible(strian::T,ttrian::RefinedTriangulation) where {T <: Triangulation}
  return is_change_possible(ttrian,strian)
end

function Geometry.best_target(strian::RefinedTriangulation,ttrian::RefinedTriangulation)
  @check is_change_possible(strian,ttrian)
  smodel = get_refined_model(strian)
  tmodel = get_refined_model(ttrian)
  a = get_model(smodel) === get_parent(tmodel)
  b = get_parent(tmodel) === smodel # tmodel = refine(smodel)
  a || b ? ttrian : strian
end

function Geometry.best_target(strian::RefinedTriangulation,ttrian::T) where {T <: Triangulation}
  @check is_change_possible(strian,ttrian)
  return strian
end

function Geometry.best_target(strian::T,ttrian::RefinedTriangulation) where {T <: Triangulation}
  @check is_change_possible(strian,ttrian)
  return ttrian
end


"""
  Given a RefinedTriangulation and a CellField defined on the parent(coarse) mesh, 
  returns an equivalent CellField on the fine mesh.
"""
function change_domain_c2f(f_coarse, ftrian::RefinedTriangulation{Dc,Dp}) where {Dc,Dp}
  model  = get_refined_model(ftrian)
  glue   = get_glue(model)
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

function CellData.change_domain(a::CellField,ttrian::RefinedTriangulation)
  strian = get_triangulation(a)
  if strian === ttrian
    return a
  end
  @assert is_change_possible(strian,ttrian)
  change_domain_c2f(a,ttrian)
end
