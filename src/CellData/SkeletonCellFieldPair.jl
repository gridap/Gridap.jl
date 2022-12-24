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
struct SkeletonCellFieldPair{
  P<:CellFieldAt, M<:CellFieldAt, T<:Triangulation} <: CellField

  cf_plus::P
  cf_minus::M
  trian::T

  function SkeletonCellFieldPair(
    cf_plus_plus::CellFieldAt, cf_minus_minus::CellFieldAt)

    @check DomainStyle(cf_plus_plus) == DomainStyle(cf_minus_minus)
    cf_plus_trian = get_triangulation(cf_plus_plus)
    cf_minus_trian =  get_triangulation(cf_minus_minus)
    @notimplementedif !(cf_plus_trian ===  cf_minus_trian)

    P = typeof(cf_plus_plus)
    M = typeof(cf_minus_minus)
    T = typeof(cf_plus_trian)
    new{P,M,T}(cf_plus_plus, cf_minus_minus, cf_plus_trian)
  end

  """
  Copy constructor
  """
  function SkeletonCellFieldPair(a::SkeletonCellFieldPair{P,M,T}) where {P,M,T}
      new{P,M,T}(copy(a.cf_plus), copy(a.cf_minus), a.trian)
  end
end

function SkeletonCellFieldPair(cf_plus::CellField, cf_minus::CellField)
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

Base.copy(a::SkeletonCellFieldPair) = SkeletonCellFieldPair(a)

function get_triangulation(a::SkeletonCellFieldPair)
  getfield(a,:trian)
end

#=
We arbitrarily choose plus side of SkeletonCellFieldPair for evaluations, as we
don't use the DomainContributions associated with Ω and Γ while performing the
sensitivities for SkeletonTriangulation (Λ) DomainContribution. This fix is for
error free evaluation of the functional when a SkeletonCellFieldPair is passed
into it. Ideally, if we could parse and extract only the Skeleton integration
terms from the functional's julia function form, this fix is not required, but
this is not trivial to do. On the positive side, since the evaluations are all
lazy and not used, this doesn't put any noticable memory or computational
overhead. Ofcourse, it is made sure that the such plus side pick doesn't happen
when the integration over the SkeletonTriangulation
=#
# If SkeletonCellFieldPair is evaluated we just pick the plus side parent
function evaluate!(cache,f::SkeletonCellFieldPair,x::CellPoint)
  _f, _x = _to_common_domain(f.plus.parent,x)
  cell_field = get_data(_f)
  cell_point = get_data(_x)
  lazy_map(evaluate,cell_field,cell_point)
end

#=
Fix for CellFieldAt{T}(parent::OperationCellField) for OperationCellField
involving args which are not FEFunctions. It creates a problem with
SkeletonCellFieldPair giving the wrong output or errors as it chooses the
SkeletonCellFieldPair directly as the parent rather the parent CellField of the
side of choice. In code terms it is the following:

CellFieldAt{T}(OperationCellField) results in calling CellFieldAt{T}(operands).
Without this override, CellFieldAt{T}(a::SkeletonCellField) results in the
CellFieldAt structure with parent as SkeletonCellField and not a.T side
CellField, which is not correct! This resulting in errors and is not the
intended behaviour to have the parent as SkeletonCellFieldPair and is not
consistent with our getproperty rules.
=#
function CellFieldAt{T}(parent::SkeletonCellFieldPair) where T
  getproperty(parent,T)
end

#=
To handle the evaluation of SkeletonCellFieldPair at BoundaryTriangulation
making it consistent with plus side choice of direct evaluation of SCFP in
general this is the case when trian of CellField is not the same as that of the
quadrature.
But we are protected from inconsistent behaviour at SkeletonTriangulation as it
fails due to similar ambiguity as the CellField at the SkeletonTrian for direct
evaluations, so get_data doesn't create any side effects.
In case of BodyFittedTriangulation there is no hit to get_data directly and
evaluate! handles it as intended
=#
get_data(a::SkeletonCellFieldPair) = get_data(a.plus)

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
