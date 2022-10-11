

struct RefinedTriangulation{Dc,Dp,T<:Triangulation} <: Triangulation{Dc,Dp}
  trian::T
end

# Wrap Triangulation API
function get_background_model(t::RefinedTriangulation)
  get_background_model(t.trian)
end

function get_grid(t::RefinedTriangulation)
  get_grid(t.trian)
end

function get_glue(t::RefinedTriangulation,::Val{d}) where d
  get_glue(t.trian,Val(d))
end

# Wrap constructors for RefinedDiscreteModel
function Triangulation(
  ::Type{ReferenceFE{d}},model::RefinedDiscreteModel,filter::AbstractArray) where d
  
  trian = Triangulation(ReferenceFE{d},model.model,filter)
  return RefinedTriangulation(trian)
end


# Domain changes
# TODO: This assumes we have the same type of triangulation on both refinement levels!
#       we might want to change this in the future when doing hybrid methods etc...

function is_change_possible(strian::RefinedTriangulation,ttrian::RefinedTriangulation)
  smodel = get_background_model(strian)
  tmodel = get_background_model(ttrian)
  a = get_model(smodel) === get_parent(tmodel) # tmodel = refine(smodel)
  b = get_parent(smodel) === get_model(tmodel) # smodel = refine(tmodel)
  return a || b
end

function is_change_possible(strian::RefinedTriangulation,ttrian::T) where {T <: Triangulation}
  smodel = get_background_model(strian)
  tmodel = get_background_model(ttrian)
  return get_parent(smodel) === tmodel # smodel = refine(tmodel)
end

function is_change_possible(strian::T,ttrian::RefinedTriangulation) where {T <: Triangulation}
  return is_change_possible(ttrian,strian)
end

function best_target(strian::RefinedTriangulation,ttrian::RefinedTriangulation)
  @check is_change_possible(strian,ttrian)
  smodel = get_background_model(strian)
  tmodel = get_background_model(ttrian)
  get_model(smodel) === get_parent(tmodel) ? ttrian : strian
end

function best_target(strian::RefinedTriangulation,ttrian::T) where {T <: Triangulation}
  @check is_change_possible(strian,ttrian)
  return strian
end
