module FESpaceConstructors

using Gridap
using Gridap.Helpers
import Gridap: FESpace


function FESpace(;kwargs...)

  reffe = _get_kwarg(:reffe,kwargs)
  @assert isa(reffe,Symbol) "For the moment, reffe can only be a symbol"


  model = _get_kwarg(:model,kwargs)
  labels = _get_kwarg(:labels,kwargs,nothing)
  conformity = _get_kwarg(:conformity,kwargs,true)

  diritags = _get_kwarg(:diritags,kwargs,Int[])
  dirimasks = _get_kwarg(:dirimasks,kwargs,nothing)

  order = _get_kwarg(:order,kwargs)

  polytope = _get_polytope(model)

  constraint = _get_kwarg(:constraint,kwargs,nothing)

  dim = celldim(model)

  fespace = nothing

  if reffe in [:Lagrangian,:QLagrangian,:PLagrangian]

    T = _get_kwarg(:valuetype,kwargs,nothing)
    if T == nothing
      error("valuetype is a mandatory keyword argument in FESpace constructor for Lagrangian reference FEs")
    end

    if _is_reffe_lagrangian_compatible_with_polytope(reffe,polytope)

      if conformity in [false, :L2]
        if labels == nothing
          fespace =  DLagrangianFESpace(T,model,order,diritags,dirimasks)
        else
          fespace =  DLagrangianFESpace(T,model,labels,order,diritags,dirimasks)
        end

      elseif conformity in [true, :default, :H1]
        if labels == nothing
          fespace = CLagrangianFESpace(T,model,order,diritags,dirimasks)
        else
          fespace = CLagrangianFESpace(T,model,labels,order,diritags,dirimasks)
        end

      else

        s = "Conformity $conformity not implemented for $reffe reference FE on a $(_poly_type(polytope))"
        error(s)

      end

    elseif reffe == :PLagrangian

      if conformity in [false, :L2]

        _reffe = PDiscRefFE(T,polytope,order)
        fespace = DiscFESpace(_reffe,model)

      else

        s = "Conformity $conformity not possible for $reffe reference FE on a $(_poly_type(polytope))"
        error(s)

      end


    else
      error("Reference element $reffe nor implemented on a $(_poly_type(polytope))")
    end

  elseif reffe == :RaviartThomas

    if ! _is_cube(polytope)
      error("RaviartThomas reference FEs can only be constructed on top of a ncube. Check your model.")
    end

    _reffe = RTRefFE(polytope,Float64,order)

    if conformity in [false, :L2]

      fespace = DiscFESpace(_reffe,model)

    elseif conformity in [true, :default, :HDiv]

      if labels == nothing
        _labels = FaceLabels(model)
      else
        _labels = labels
      end

      grid = Grid(model)
      trian = Triangulation(grid)
      graph = GridGraph(model)

      # fespace = ConformingFESpace(_reffe,trian,graph,_labels,diritags)
      fespace = DivConformingFESpace(_reffe,trian,graph,_labels,diritags)

    else

      s = "Conformity $conformity not possible for a $reffe reference FE"
      error(s)

    end

  elseif reffe == :Nedelec

    if ! _is_cube(polytope)
      error("Nedelec reference FEs can only be constructed on top of a ncube. Check your model.")
    end

    _reffe = NedelecRefFE(polytope,Float64,order)

    if conformity in [false, :L2]

      fespace = DiscFESpace(_reffe,model)

    elseif conformity in [true, :default, :HCurl]

      if labels == nothing
        _labels = FaceLabels(model)
      else
        _labels = labels
      end

      grid = Grid(model)
      trian = Triangulation(grid)
      graph = GridGraph(model)

      fespace = ConformingFESpace(_reffe,trian,graph,_labels,diritags)

    else

      s = "Conformity $conformity not possible for a $reffe reference FE"
      error(s)

    end

  else
    error("Reference element $reffe not implemented")
  end

  @assert fespace != nothing

  if constraint == nothing
    return fespace

  elseif constraint == :zeromean
    return ZeroMeanFESpace(fespace,order)

  else
    error("Unknown constraint value $constraint")

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

function _get_polytope(model)

  polys = CellPolytopes(Grid(model))

  @notimplementedif ! isa(polys,ConstantCellValue)

  polys.value

end

function _is_reffe_lagrangian_compatible_with_polytope(reffe,polytope)
  v = (reffe == :Lagrangian)
  v = v || (reffe == :QLagrangian && _is_cube(polytope))
  v = v || (reffe == :PLagrangian && _is_simplex(polytope))
  v
end

_is_cube(p) = all(p.extrusion.array .== HEX_AXIS)

_is_simplex(p) = all(p.extrusion.array .== TET_AXIS)

function _poly_type(p)

  if _is_cube(p)
    t = :ncube
  elseif _is_simplex(p)
    t = :simplex
  else
    @notimplemented
  end

  return t

end



end # module
