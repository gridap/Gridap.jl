
"""
    divergence(f)

Abstract divergence operator, formally equivalent to `f -> ∇⋅f`.
"""
divergence(f) = Operation(tr)(∇(f))

"""
    DIV(f)

Reference space divergence.
"""
function DIV(f)
  @notimplemented "DIV operator is only defined for certain type of operands"
end

function symmetric_gradient end

"""
    symmetric_gradient(f)
    ε(f)

Abstract symmetric gradient operator, formally equivalent to `f -> ½(∇f + (∇f)ᵀ)`.
"""
symmetric_gradient(f) = Operation(symmetric_part)(gradient(f))
const ε = symmetric_gradient

"""
    skew_symmetric_gradient(f)

Abstract skew symmetric gradient operator, formally equivalent to `f -> ½(∇f - (∇f)ᵀ)`.
"""
skew_symmetric_gradient(f) = Operation(skew_symmetric_part)(gradient(f))

"""
    curl(f)

Abstract curl operator, formally equivalent to
- `f -> ∂₁f₂ - ∂₂f₁` for 2D vector functions, or
- `f -> ∇×f` for 3D vector functions.
"""
curl(f) = Operation(grad2curl)(∇(f))

"""
    grad2curl(∇f)

Return
- `∇f[1,2] - ∇f[2,1]` for 2×2 input tensor, or
- `VectorValue(∇f[2,3] - ∇f[3,2], ∇f[3,1] - ∇f[1,3], ∇f[1,2] - ∇f[2,1])` for 3×3 input tensor.
"""
function grad2curl(∇u::TensorValue{2,2})
  ∇u[1,2] - ∇u[2,1]
end

function grad2curl(∇u::TensorValue{3,3})
  c1 = ∇u[2,3] - ∇u[3,2]
  c2 = ∇u[3,1] - ∇u[1,3]
  c3 = ∇u[1,2] - ∇u[2,1]
  VectorValue(c1,c2,c3)
end

function laplacian end
const Δ = laplacian

"""
    laplacian(f)
    Δ(f)

Abstract laplacian operator, equivalent to `tr(∇∇(f))`.
"""
function laplacian(f)
  # g = gradient(f)
  # divergence(g)
  tr(∇∇(f))
end

"""
    dot(∇,f)
    ∇⋅f

Equivalent to `divergence(f)`.
"""
dot(::typeof(∇),f::Field) = divergence(f)
dot(::typeof(∇),f::Function) = divergence(f)

function (*)(::typeof(∇),f::Field)
  msg = "Syntax ∇*f has been removed, use ∇⋅f (\\nabla \\cdot f) instead"
  error(msg)
end

function (*)(::typeof(∇),f::Function)
  msg = "Syntax ∇*f has been removed, use ∇⋅f (\\nabla \\cdot f) instead"
  error(msg)
end

"""
    outer(∇,f)
    ∇⊗f

Equivalent to `gradient(f)`.
"""
outer(::typeof(∇),f::Field) = gradient(f)
outer(::typeof(∇),f::Function) = gradient(f)

"""
    outer(f,∇)
    f⊗∇

Equivalent to `transpose(gradient(f))`.
"""
outer(f::Field,::typeof(∇)) = transpose(gradient(f))
outer(f::Function,::typeof(∇)) = transpose(gradient(f))

"""
    cross(∇,f)
    ∇×f

Equivalent to `curl(f)`.
"""
cross(::typeof(∇),f::Field) = curl(f)
cross(::typeof(∇),f::Function) = curl(f)

_extract_grad_diag(x::TensorValue) = diag(x)
_extract_grad_diag(x) = @notimplemented

function Base.broadcasted(::typeof(*),::typeof(∇),f)
  g = ∇(f)
  Operation(_extract_grad_diag)(g)
end

function Base.broadcasted(::typeof(*),::typeof(∇),f::Function)
  Base.broadcasted(*,∇,GenericField(f))
end

# Shifted ∇

struct ShiftedNabla{N,T}
  v::VectorValue{N,T}
end

(+)(::typeof(∇),v::VectorValue) = ShiftedNabla(v)
(+)(v::VectorValue,::typeof(∇)) = ShiftedNabla(v)
(-)(::typeof(∇),v::VectorValue) = ShiftedNabla(-v)

function (s::ShiftedNabla)(f)
  Operation((a,b)->a+s.v⊗b)(gradient(f),f)
end

(s::ShiftedNabla)(f::Function) = s(GenericField(f))

function return_value(k::Broadcasting{<:ShiftedNabla},f)
  s = k.f
  g = Broadcasting(∇)(f)
  Broadcasting(Operation((a,b)->a+s.v⊗b))(g,f)
end

function evaluate!(cache,k::Broadcasting{<:ShiftedNabla},f)
  s = k.f
  g = Broadcasting(∇)(f)
  Broadcasting(Operation((a,b)->a+s.v⊗b))(g,f)
end

dot(s::ShiftedNabla,f) = Operation(tr)(s(f))
outer(s::ShiftedNabla,f) = s(f)
outer(f,s::ShiftedNabla) = transpose(s(f))
cross(s::ShiftedNabla,f) = Operation(grad2curl)(s(f))

dot(s::ShiftedNabla,f::Function) = dot(s,GenericField(f))
outer(s::ShiftedNabla,f::Function) = outer(s,GenericField(f))
outer(f::Function,s::ShiftedNabla) = outer(GenericField(f),s)
cross(s::ShiftedNabla,f::Function) = cross(s,GenericField(f))

function Base.broadcasted(::typeof(*),s::ShiftedNabla,f)
  g = s(f)
  Operation(_extract_grad_diag)(g)
end

function Base.broadcasted(::typeof(*),s::ShiftedNabla,f::Function)
  Base.broadcasted(*,s,GenericField(f))
end

# Variance consistent differential operators, a.k.a. Differential geometry operators

"""
    exterior_derivative(ω)
    𝑑(ω)

Exterior derivative of a differential K-form `dω`, is a (K+1)-form.
"""
exterior_derivative(ω) = Operation(partials_to_exteriorderivative)(partial_derivatives(ω))
const 𝑑 = exterior_derivative

"""
    ∧(∇,f)
    ∇∧f

Equivalent to `exterior_derivative(f)`.
"""
∧(::typeof(∇),f::Field) = exterior_derivative(f)
∧(::typeof(∇),f::Function) = exterior_derivative(f)

"""
    lie_derivative(v, ω)
    𝓛(v,ω)

Lie derivative of the the field `ω` with respect to the vector field `v`. If
`ω` is a k-form, 𝓛_v ω = d(ι_v ω) + ι_v(dω) by Cartan's magic formula.

`v` is a contravariant vector field and `ω` a scalar or tensor field with
variance of any order. `𝓛_v ω` is of the same type as `ω`, modulo component type
promotion with that of `v`.
"""
lie_derivative(v,ω) = Operation(partials_to_liederivative)(v,ω,partial_derivatives(v),partial_derivatives(ω))
const 𝓛 = lie_derivative

"""
    covariant_derivative(a::CellField)
    covariant_derivative(v::CellField, a::CellField)
    covariant_derivative([v::F, ]a::F, g)
    cov∇([v, ]a[, g])

Covariant derivative of the field `a` for the Levi-Civita connection of the
metric `g`, along the (contravariant) vector field `v` if given. `cov∇(v, a)` is
of the same type as `a`, `cov∇(a)` has the type of `∂(a)`.

It depends on the metric, which is automatically fed to `a` if it is a
`CellField`, or needed explicitly for `Field` or a `Function` arguments. `g`
may be `LinearAlgebra.I`. `cov∇(a, I)` is the same as `∂`.
"""
function covariant_derivative end

covariant_derivative(a, ::UniformScaling) = partial_derivatives(a)

covariant_derivative(v, a, ::UniformScaling) =
  Operation(partials_to_flat_covariantderivative)(v, partial_derivatives(a))

const cov∇ = covariant_derivative

"""
    codifferential(ω::CellField)
    codifferential(ω::F, g)
    δ(ω[, g])

Codifferential of a differential K-form `ω`, is a (K-1)-form:
  ``δω = s(-1)^{D(K-1)+1} ⋆ d ⋆ ω``
with ``s = sgn(det(g))``, ``+1`` for Riemanian metric, ``-1`` for Lorentzian metric.

It depends on the metric, which is automatically fed to `ω` if it is a cell
field. For a `Function` or `Field` `ω::F`, `g` must be provided, but it can be
`LinearAlgebra.I` (Euclidean metric).
"""
function codifferential end

codifferential(ω, ::UniformScaling) =
  Operation(partials_to_flat_codifferential)(partial_derivatives(ω))

function codifferential(ω, g)
  Γ = Operation(partials_to_Christoffel_2)(∂(g), g)
  Operation(partials_to_codifferential)(∂(ω), ω, Γ, g)
end

const δ = codifferential

"""
    bochner_laplacian(a::F, g)
    ΔB(a, g)

Bochner laplacian of the `Field` or `Function` `a::F`, `ΔB = ∇*∇ = -tr ∇∇`,
positive semi-definite, so minus the connection laplacian. `ΔB a` is of the same
type as `a`, and is minus the Laplace-Beltrami operator [`laplacian`](@ref) on a
scalar field.

The metric `g` is required, it may be of the type of `a` or `LinearAlgebra.I`.
"""
function bochner_laplacian end

bochner_laplacian(a, ::UniformScaling) =
  Operation(second_partials_to_flat_bochnerlaplacian)(∂∂(a))

bochner_laplacian(a, g) =
  Operation(second_partials_to_bochnerlaplacian)(∂∂(a), ∂(a), a, ∂∂(g), ∂(g), g)

const ΔB = bochner_laplacian

"""
    hodge_laplacian(ω::F, g)
    ΔdR(ω, g)

Hodge laplacian of a differential K-form valued `Field` or `Function` `ω`, `ΔdR
= dδ + δd`, is a K-form.

The metric `g` is required, it may be of the type of `a` or `LinearAlgebra.I`.
"""
function hodge_laplacian end

hodge_laplacian(ω, ::UniformScaling) =
  Operation(second_partials_to_flat_hodgelaplacian)(∂∂(ω))

hodge_laplacian(ω, g) =
  Operation(second_partials_to_hodgelaplacian)(∂∂(ω), ∂(ω), ω, ∂∂(g), ∂(g), g)

const ΔdR = hodge_laplacian

"""
    koszul(ω)

Koszul differential of a differential K-form `ω`, is a (K-1)-form.

`κ_x(ω) = ι_x(ω)` contracts `ω` with the evaluation point, which the identity
field supplies as a second argument.
"""
koszul(ω::Field) = Operation(koszul)(GenericField(identity),ω)
koszul(ω::Function) = koszul(GenericField(ω))

"""
    d_0form(f) = to_1form(∇(f))

Discrete exterior derivative of a scalar (0-form) cell field.
Returns a `CoVectorValue{D}`-valued `OperationCellField`.
"""
d_0form(f) = to_1form(∇(f))

"""
    d_1form(f) = Operation(grad_to_2form)(∇(f))

Discrete exterior derivative of a `VectorValue`'d (1-form) CellField
(e.g., from a Nédélec FESpace).
Returns a `ExteriorFormValue{2,D}`-valued `OperationCellField`.
"""
d_1form(f) = Operation(grad_to_2form)(∇(f))

# Broadcasting of differential operators using gradient
for (diff_op, op) in (
  (:divergence, :tr), (:symmetric_gradient,:symmetric_part),
  (:skew_symmetric_gradient,:skew_symmetric_part),
  (:curl,:grad2curl), (:d_0form,:to_1form), (:d_1form,:grad_to_2form),
)
  @eval begin
    function return_value(::Broadcasting{typeof($diff_op)},f)
      Broadcasting(Operation($op))(Broadcasting(∇)(f))
    end
    function evaluate!(cache,::Broadcasting{typeof($diff_op)},f)
      Broadcasting(Operation($op))(Broadcasting(∇)(f))
    end
  end
end

# Broadcasting of differential operators using partial_derivatives.

function return_value(::Broadcasting{typeof(exterior_derivative)},f)
  Broadcasting(Operation(partials_to_exteriorderivative))(Broadcasting(partial_derivatives)(f))
end

function evaluate!(cache,::Broadcasting{typeof(exterior_derivative)},f)
  Broadcasting(Operation(partials_to_exteriorderivative))(Broadcasting(partial_derivatives)(f))
end

# `δ` reads the metric, which a field array carries no more than a field does,
# so its broadcast takes the metric argument the operator itself takes.
return_value(k::Broadcasting{typeof(codifferential)},f,metric) =
  evaluate!(nothing,k,f,metric)

function evaluate!(cache,::Broadcasting{typeof(codifferential)},f,::UniformScaling)
  Broadcasting(Operation(partials_to_flat_codifferential))(Broadcasting(∂)(f))
end

function evaluate!(cache,::Broadcasting{typeof(codifferential)},f,g)
  Γ = Broadcasting(Operation(partials_to_Christoffel_2))(Broadcasting(∂)(g), g)
  Broadcasting(Operation(partials_to_codifferential))(Broadcasting(∂)(f), f, Γ, g)
end

for (op, flat_kernel, kernel) in (
  (:bochner_laplacian, :second_partials_to_flat_bochnerlaplacian, :second_partials_to_bochnerlaplacian),
  (:hodge_laplacian,   :second_partials_to_flat_hodgelaplacian,   :second_partials_to_hodgelaplacian),
)
  @eval begin
    return_value(k::Broadcasting{typeof($op)}, f, metric) =
      evaluate!(nothing, k, f, metric)

    function evaluate!(cache,::Broadcasting{typeof($op)}, f, ::UniformScaling)
      Broadcasting(Operation($flat_kernel))(Broadcasting(∂∂)(f))
    end

    function evaluate!(cache,::Broadcasting{typeof($op)}, f, g)
      Broadcasting(Operation($kernel))(
        Broadcasting(∂∂)(f), Broadcasting(∂)(f), f,
        Broadcasting(∂∂)(g), Broadcasting(∂)(g), g)
    end
  end
end

return_value(k::Broadcasting{typeof(covariant_derivative)}, a, args...) =
  evaluate!(nothing, k, a, args...)

evaluate!(cache, ::Broadcasting{typeof(covariant_derivative)}, a, ::UniformScaling) =
  Broadcasting(partial_derivatives)(a)

function evaluate!(cache, ::Broadcasting{typeof(covariant_derivative)}, v, a, ::UniformScaling)
  Broadcasting(Operation(partials_to_flat_covariantderivative))(
    v, Broadcasting(partial_derivatives)(a))
end

