# TODO Move the content of this file into src/Exports.jl

macro publish_gridapodes(name)
  quote
    using Gridap.ODEs: $name; export $name
  end
end

@publish_gridapodes ∂t
@publish_gridapodes ∂tt
# @publish_gridapodes BackwardEuler
# @publish_gridapodes ForwardEuler
# @publish_gridapodes ThetaMethod
# @publish_gridapodes MidPoint
# @publish_gridapodes RungeKutta
# @publish_gridapodes IMEXRungeKutta
# @publish_gridapodes Newmark
# @publish_gridapodes GeneralizedAlpha

@publish_gridapodes TransientTrialFESpace
@publish_gridapodes TransientMultiFieldFESpace
@publish_gridapodes TransientFEOperator
# @publish_gridapodes TransientMassLinearFEOperator
