"""
  SkeletonCellFieldPair is a special construct for allowing uh.plus
  and uh.minus to be two different CellFields. In particular, it is
  useful when we need one of the CellFields to be the dualized
  version of the other for performing ForwardDiff AD of a Skeleton
  integration DomainContribution wrt to the degrees of freedom
  of the CellField, plus and minus sensitivities done separately,
  so as to restrict the interaction between the dual numbers.

  It takes in two CellFields and stores plus version of CellFieldAt
  of the first CellField and minus version of CellFieldAt of the
  second the CellField. SkeletonCellFieldPair is associated with
  same triangulation as that of the CellFields (we check if the
  triangulations of both CellFields match)

  SkeletonCellFieldPair is an internal convenience artifact/construct
  to aid in dualizing plus and minus side around a Skeleton face separately
  to perform the sensitivity of degrees of freedom of cells sharing the
  Skeleton face, without the interaction dual numbers of the two cells.
  The user doesn't have to deal with this construct anywhere when
  performing AD of functionals involving integration over Skeleton faces
  using the public API.
"""
struct SkeletonCellFieldPair{P<:CellFieldAt, M<:CellFieldAt} <: CellField
  cf_plus::P
  cf_minus::M

  function SkeletonCellFieldPair(cf_plus_plus,cf_minus_minus)
    P = typeof(cf_plus_plus)
    M = typeof(cf_minus_minus)
    new{P,M}(cf_plus_plus, cf_minus_minus)
  end
end

function SkeletonCellFieldPair(
  cf_plus::CellField,
  cf_minus::CellField,
  strian::SkeletonTriangulation)

  @check DomainStyle(cf_plus) == DomainStyle(cf_minus)
  cf_plus_trian = get_triangulation(cf_plus)
  cf_minus_trian =  get_triangulation(cf_minus)
  @notimplementedif !(cf_plus_trian ===  cf_minus_trian)

  @check num_dims(cf_plus_trian) == num_dims(strian) + 1
  @check get_background_model(cf_plus_trian) === get_background_model(strian)

  cf_plus_plus  = cf_plus.plus
  cf_minus_minus = cf_minus.minus

  SkeletonCellFieldPair(cf_plus_plus,cf_minus_minus)
end

function Base.getproperty(a::SkeletonCellFieldPair,sym::Symbol)
  if sym in (:⁺,:plus)
    getfield(a,:cf_plus)
  elseif sym in (:⁻, :minus)
    getfield(a,:cf_minus)
  else
    getfield(a,sym)
  end
end

function DomainStyle(a::SkeletonCellFieldPair)
  DomainStyle(getfield(a,:cf_plus))
end

function get_triangulation(a::SkeletonCellFieldPair)
  get_triangulation(getfield(a,:cf_plus))
end

#=
The below change_domain has been overloaded for SkeletonCellFieldPair
to work without errors for integrations over other triangulation - Ω
(BodyFittedTriangulation) and Γ (BoundaryTriangulation).

We arbitrarily choose plus side of SkeletonCellFieldPair for evaluations, as
we don't use the DomainContributions associated with Ω and Γ while performing
the sensitivities for SkeletonTriangulation (Λ) DomainContribution. This
fix is for error free evaluation of the functional when a SkeletonCellFieldPair
is passed into it. Ideally, if we could parse and extract only the Skeleton
integration terms from the functional's julia function form, this fix is not
required, but this is not trivial to do. On the positive side, since the
evaluations are all lazy and not used, this doesn't put any noticable memory
or computational overhead.
=#
function change_domain(a::SkeletonCellFieldPair,target_trian::Triangulation,target_domain::DomainStyle)
  change_domain(a.plus,DomainStyle(a),target_trian,target_domain)
end

function Base.propertynames(a::SkeletonCellFieldPair, private::Bool=false)
  (fieldnames(typeof(a))...,:⁺,:plus,:⁻,:minus)
end

function Base.show(io::IO,::MIME"text/plain",f::SkeletonCellFieldPair)
  show(io,f)
end

# derivatives of the SkeletonCellFieldPair

function gradient(a::SkeletonCellFieldPair)
  grad_cf_plus_plus = gradient(getfield(a,:cf_plus))
  grad_cf_minus_minus = gradient(getfield(a,:cf_minus))
  SkeletonCellFieldPair(grad_cf_plus_plus,grad_cf_minus_minus)
end

function ∇∇(a::SkeletonCellFieldPair)
  hess_cf_plus_plus = ∇∇(getfield(a,:cf_plus))
  hess_cf_minus_minus = ∇∇(getfield(a,:cf_minus))
  SkeletonCellFieldPair(hess_cf_plus_plus,hess_cf_minus_minus)
end
