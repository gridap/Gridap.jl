
# Quadrature rule

"""
"""
struct CellQuadrature{DDS,IDS} <: CellDatum
  cell_quad::AbstractArray{<:Quadrature}
  cell_point::AbstractArray{<:AbstractArray{<:Point}}
  cell_weight::AbstractArray{<:AbstractArray{<:Real}}
  trian::Triangulation
  data_domain_style::DDS
  integration_domain_style::IDS
end

# Old constructor (for backward compatibility)
function CellQuadrature(
  cell_quad::AbstractArray{<:Quadrature},
  cell_point::AbstractArray{<:AbstractArray{<:Point}},
  cell_weight::AbstractArray{<:AbstractArray{<:Real}},
  trian::Triangulation,
  data_domain_style::DomainStyle)
  CellQuadrature(cell_quad,cell_point,cell_weight,trian,data_domain_style,PhysicalDomain())
end

function CellQuadrature(trian::Triangulation,quad::Tuple{<:QuadratureName,Any,Any};
                        data_domain_style::DomainStyle=ReferenceDomain(),
                        integration_domain_style::DomainStyle=PhysicalDomain())
  name, args, kwargs = quad
  cell_quad = Quadrature(trian,name,args...;kwargs...)
  CellQuadrature(trian,cell_quad,data_domain_style,integration_domain_style)
end

function CellQuadrature(trian::Triangulation,degree::Integer;
                        data_domain_style::DomainStyle=ReferenceDomain(),
                        integration_domain_style::DomainStyle=PhysicalDomain(),kwargs...)
  cell_quad = Quadrature(trian,degree;kwargs...)
  CellQuadrature(trian,cell_quad,data_domain_style,integration_domain_style)
end

function CellQuadrature(trian::Triangulation,quad::Quadrature;
                        data_domain_style::DomainStyle=ReferenceDomain(),
                        integration_domain_style::DomainStyle=PhysicalDomain())
  cell_quad = expand_cell_data([quad],get_cell_type(trian))
  CellQuadrature(trian,cell_quad,data_domain_style,integration_domain_style)
end

function CellQuadrature(trian::Triangulation,cell_quad::AbstractVector{<:Quadrature};
                        data_domain_style::DomainStyle=ReferenceDomain(),
                        integration_domain_style::DomainStyle=PhysicalDomain())
  CellQuadrature(trian,cell_quad,data_domain_style,integration_domain_style)
end

function CellQuadrature(trian::Triangulation,
                        cell_quad::AbstractVector{<:Quadrature},
                        dds::DomainStyle,ids::DomainStyle)
  cell_point = lazy_map(get_coordinates,cell_quad)
  cell_weight = lazy_map(get_weights,cell_quad)
  CellQuadrature(cell_quad,cell_point,cell_weight,trian,dds,ids)
end

function CellQuadrature(trian::AppendedTriangulation,degree1,degree2;
                        data_domain_style::DomainStyle=ReferenceDomain(),
                        integration_domain_style::DomainStyle=PhysicalDomain(),kwargs...)
  quad1 = CellQuadrature(trian.a,degree1;
                         data_domain_style=data_domain_style,
                         integration_domain_style=integration_domain_style,kwargs...)
  quad2 = CellQuadrature(trian.b,degree2;
                         data_domain_style=data_domain_style,
                         integration_domain_style=integration_domain_style,kwargs...)
  lazy_append(quad1,quad2,trian)
end

function CellQuadrature(trian::AppendedTriangulation,degree::Integer;
                        data_domain_style::DomainStyle=ReferenceDomain(),
                        integration_domain_style::DomainStyle=PhysicalDomain(),kwargs...)
  CellQuadrature(trian,degree,degree;
                 data_domain_style=data_domain_style,
                 integration_domain_style=integration_domain_style,kwargs...)
end

function CellQuadrature(trian::AppendedTriangulation,
                        quad::Tuple{<:QuadratureName,Any,Any};
                        data_domain_style::DomainStyle=ReferenceDomain(),
                        integration_domain_style::DomainStyle=PhysicalDomain())
  CellQuadrature(trian,quad,quad;
                 data_domain_style=data_domain_style,
                 integration_domain_style=integration_domain_style)
end

function CellQuadrature(trian::AppendedTriangulation,quad::Quadrature;
                        data_domain_style::DomainStyle=ReferenceDomain(),
                        integration_domain_style::DomainStyle=PhysicalDomain())
  CellQuadrature(trian,quad,quad;
                 data_domain_style=data_domain_style,
                 integration_domain_style=integration_domain_style)
end

@deprecate(
  CellQuadrature(trian::Triangulation,degree,ids::DomainStyle),
  CellQuadrature(trian::Triangulation,degree;data_domain_style=ReferenceDomain(),integration_domain_style=ids)
)

@deprecate(
  CellQuadrature(trian::AppendedTriangulation,degree1,degree2,ids::DomainStyle),
  CellQuadrature(trian::AppendedTriangulation,degree1,degree2;data_domain_style=ReferenceDomain(),integration_domain_style=ids)
)

function lazy_append(
  quad1::CellQuadrature,
  quad2::CellQuadrature,
  trian::AppendedTriangulation=lazy_append(quad1.trian,quad2.trian))

  @notimplementedif quad1.data_domain_style != quad2.data_domain_style
  @notimplementedif quad1.integration_domain_style != quad2.integration_domain_style
  cell_quad = lazy_append(quad1.cell_quad,quad2.cell_quad)
  cell_point = lazy_append(quad1.cell_point,quad2.cell_point)
  cell_weight = lazy_append(quad1.cell_weight,quad2.cell_weight)
  CellQuadrature(cell_quad,cell_point,cell_weight,trian,
                 quad1.data_domain_style,quad1.integration_domain_style)
end

get_data(f::CellQuadrature) = f.cell_quad
get_triangulation(f::CellQuadrature) = f.trian
DomainStyle(::Type{CellQuadrature{DDS,IDS}}) where {DDS,IDS} = DDS()

function change_domain(a::CellQuadrature,::ReferenceDomain,::PhysicalDomain)
  @notimplemented
end

function change_domain(a::CellQuadrature,::PhysicalDomain,::ReferenceDomain)
  @notimplemented
end

function get_cell_points(a::CellQuadrature)
  CellPoint(a.cell_point,a.trian,a.data_domain_style)
end

function integrate(f::CellField, quad::CellQuadrature)
  trian_f = get_triangulation(f)
  trian_x = get_triangulation(quad)

  msg = """\n
    Your are trying to integrate a CellField using a CellQuadrature defined on incompatible
    triangulations. Verify that either the two objects are defined in the same triangulation
    or that the triangulaiton of the CellField is the background triangulation of the CellQuadrature.
    """
  @check is_change_possible(trian_f,trian_x) msg

  b = change_domain(f, quad.trian, quad.data_domain_style)
  x = get_cell_points(quad)

  T = eltype(testitem(b(x)))
  bx = _integrand_coefficient(T, b, quad.trian)(x)

  if quad.data_domain_style == quad.integration_domain_style
    lazy_map(IntegrationMap(), bx, quad.cell_weight)
  elseif quad.data_domain_style == ReferenceDomain() &&
            quad.integration_domain_style == PhysicalDomain()
    cell_Jt = get_cell_jac(quad.trian)
    cell_Jtx = lazy_map(evaluate, cell_Jt, quad.cell_point)
    # The contribution of the triangulation's intrinsic metric is added to the
    # weights rather than the integrand, which is more efficient for cell-wise
    # constant metric.
    vol_weights = _metric_weights(T, x, quad, get_metric_data(quad.trian),)
    lazy_map(IntegrationMap(), bx, vol_weights, cell_Jtx)
  else
    @notimplemented
  end
end

_integrand_coefficient(::Type, b::CellField, trian::Triangulation) = b
function _integrand_coefficient(
  ::Type{<:ExteriorFormValue{K,Dw}}, b::CellField, trian::Triangulation) where {K,Dw}

  # 0-form: same as integral of scalar function
  K == 0 && return from_0form(b)

  # b = b_I dxᴵ where length(I) = manifold dimension:
  # ∫( b ) = ∫( b_I )dΩ (the latter being the usual integral).
  K == Dw == num_cell_dims(trian) && return flat_vol_coeff(b)
  # the "vol_coeff" name here is misleading: this is just the form components
  # in the standard basis, the volume form contribution is delt with in
  # _metric_weights, see [`flat_vol_coeff`](@ref).

  @unreachable """\n
  A $K-form of $(Int(Dw)) dimensions has no integral over a $(num_cell_dims(trian))-dimensional
  domain. The form valued fields that can be integrated are:
  - a 0-form, which is integrated like a scalar function,
  - a top form of the domain, whose single coefficient is what the integral of a form means, have one.
  Look at using an inner product such as `alt_inner(u,v) === u ∧ ⋆v` or
  `inner`, or integrate its components as those of the vector `from_Kform(u)`.
  """
end

# The volume element of the integral of a function over a triangulation with
# metric `g`, pulled back in reference space, is `meas(Jt)·√|det g|` where `Jt`
# is the cell_map Jacobian from the triangulation's element to the reference
# element.  meas(Jt) is added by the integration map, √|det g| is added here.
# For a volume form, √|det g| is accounted for naturally so must not be added.
function _metric_weights(::Type{T}, x::CellPoint, quad::CellQuadrature, metric_data) where T
  if (metric_data[1] isa UniformScaling) || (T <: ExteriorFormValue && !iszero(get_form_degree(T)))
    return quad.cell_weight
  end
  # _integrand_coefficient guaranties this is a scalar integrand
  sqrt_abs_g = metric_data[3]
  lazy_map(Broadcasting(*), quad.cell_weight, sqrt_abs_g(x))
end

"""
    integrate(integrand, dΩ::Measure)
    integrate(integrand, dΩ::CellQuadrature)
    ∫(integrand) * dΩ
    ∫(integrand)dΩ

High level integral definition API. The `integrand` can be a `CellField` such
as an `FEFunction`, or something that the [`CellField`](@ref) factory
accepts as a first argument like a closed form julia `Function` or constant
`Number`.

The definition of the integral depends on the return type of the integrand. If
it is a scalar valued field, it is the standard (Riemannian, Lebesgue)
integral. If the measure's triangulation declares a non-identity metric `g`
tensor via [`get_cell_metric`](@ref), the volume element √|det `g`| is
automatically added to the integrand.

If the value type is a tensor type different from a differential form value,
the integral is performed component-wise.

If it is a differential form, that is its value type is `ExteriorFormValue{K,D}`, then:
- if `K=0`, the integral is defined as above,
- if `K=D` and `D` is the manifold dimension of the triangulation, proper
  differential form integral is used: the single form coefficient is integrated
  without adding the √|det `g`| term, since the form component already contains it,
- else, the integral is not defined.

As a consequence, integrating a `0`- or `D`-form `ω` can equivalently be done
using `∫( ω )dΩ` or `∫( ⋆(ω) )dΩ`, since `⋆` acts on `0`- and `D`-forms by
multiplying by √|det `g`| and √|det `g`⁻¹|, respectively. See also
[`alt_inner`](@ref).

The syntaxic sugar `∫(integrand)dΩ` relies on wrapping the integrand into an
[`Integrand`](@ref Integrand).
"""
function integrate(a,quad::CellQuadrature)
  b = CellField(a,quad.trian,quad.data_domain_style)
  integrate(b,quad)
end

# Some syntactic sugar

"""
    struct Integrand object end
    ∫(object)

Generic placeholder for a quantity `object` to be integrated against a [`CellQuadrature`](@ref).
"""
struct Integrand
  object
end

const ∫ = Integrand

(*)(a::Integrand,b::CellQuadrature) = integrate(a.object,b)
(*)(b::CellQuadrature,a::Integrand) = integrate(a.object,b)

# Cell measure

"""
    get_cell_measure(trian)
    get_cell_measure(strian, ttrian)

where all arguments are [`Triangulation`](@ref)s, returns a vector containing
the volume of each cell of `[t]trian`.
"""
function get_cell_measure(trian::Triangulation)
  quad = CellQuadrature(trian,0)
  cell_to_dV = integrate(1,quad)
end

function get_cell_measure(strian::Triangulation,ttrian::Triangulation)
  scell_measure = get_cell_measure(strian)
  move_contributions(scell_measure,strian,ttrian) |> collect
end

#"""
#Contributions added to the cells of the background Triangulation.
#"""
#function get_cell_measure(trian::Triangulation)
#  quad = CellQuadrature(trian,0)
#  cell_to_dV = integrate(1,quad)
#  model = get_background_model(trian)
#  bgtrian = Triangulation(model)
#  D = num_cell_dims(model)
#  glue = get_glue(trian,Val(D))
#  cell_to_bgcell = glue.tface_to_mface
#  bgcell_to_dV = zeros(num_cells(bgtrian))
#  _meas_K_fill!(bgcell_to_dV,cell_to_dV,cell_to_bgcell)
#  bgcell_to_dV
#end
#
#function _meas_K_fill!(bgcell_to_dV,cell_to_dV,cell_to_bgcell)
#  cache = array_cache(cell_to_dV)
#  for (cell, bgcell) in enumerate(cell_to_bgcell)
#    dV = getindex!(cache,cell_to_dV,cell)
#    bgcell_to_dV[bgcell] += dV
#  end
#end
