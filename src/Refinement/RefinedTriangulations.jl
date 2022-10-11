

struct RefinedTriangulation{Dc,Dp,T<:Triangulation} <: Triangulation{Dc,Dp}
  trian::T

  function RefinedTriangulation(trian::Triangulation{Dc,Dp}) where {Dc,Dp}
    T = typeof(trian)
    return new{Dc,Dp,T}(trian)
  end
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
  v = view(t.trian)
  return RefinedTriangulation(v)
end

# Wrap constructors for RefinedDiscreteModel
function Geometry.Triangulation(
  ::Type{ReferenceFE{d}},model::RefinedDiscreteModel,filter::AbstractArray) where d
  
  trian = Triangulation(ReferenceFE{d},model.model,filter)
  return RefinedTriangulation(trian)
end

function Geometry.Triangulation(
  ::Type{ReferenceFE{d}},model::RefinedDiscreteModel,labels::FaceLabeling;kwargs...) where d
  trian = Triangulation(ReferenceFE{d},get_model(model),labels;kwargs...)
  return RefinedTriangulation(trian)
end

function Geometry.Triangulation(trian::RefinedTriangulation,args...;kwargs...)
  t = Triangulation(trian.trian)
  return RefinedTriangulation(t)
end

# Domain changes
# TODO: This assumes we have the same type of triangulation on both refinement levels!
#       we might want to change this in the future when doing hybrid methods etc...

function Geometry.is_change_possible(strian::RefinedTriangulation,ttrian::RefinedTriangulation)
  smodel = get_background_model(strian)
  tmodel = get_background_model(ttrian)
  a = get_model(smodel) === get_parent(tmodel) # tmodel = refine(smodel)
  b = get_parent(smodel) === get_model(tmodel) # smodel = refine(tmodel)
  return a || b
end

function Geometry.is_change_possible(strian::RefinedTriangulation,ttrian::T) where {T <: Triangulation}
  smodel = get_background_model(strian)
  tmodel = get_background_model(ttrian)
  return get_parent(smodel) === tmodel # smodel = refine(tmodel)
end

function Geometry.is_change_possible(strian::T,ttrian::RefinedTriangulation) where {T <: Triangulation}
  return is_change_possible(ttrian,strian)
end

function Geometry.best_target(strian::RefinedTriangulation,ttrian::RefinedTriangulation)
  @check is_change_possible(strian,ttrian)
  smodel = get_background_model(strian)
  tmodel = get_background_model(ttrian)
  get_model(smodel) === get_parent(tmodel) ? ttrian : strian
end

function Geometry.best_target(strian::RefinedTriangulation,ttrian::T) where {T <: Triangulation}
  @check is_change_possible(strian,ttrian)
  return strian
end
