
"""
    divergence(f)

Abstract divergence operator, formally equivalent to `f -> ∇⋅f`.
"""
divergence(f) = Operation(tr)(∇(f))

function evaluate!(cache,::Broadcasting{typeof(divergence)},f)
  Broadcasting(Operation(tr))(Broadcasting(∇)(f))
end

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

Abstract symmetric gradient operator, formally equivalent to `f -> ½(∇f + (∇f)ᵀ)`.
"""
symmetric_gradient(f) = Operation(symmetric_part)(gradient(f))

function evaluate!(cache,::Broadcasting{typeof(symmetric_gradient)},f)
  Broadcasting(Operation(symmetric_part))(Broadcasting(∇)(f))
end

"""
    const ε = symmetric_gradient

Alias for the [`symmetric_gradient`](@ref).
"""
const ε = symmetric_gradient

"""
    skew_symmetric_gradient(f)

Abstract skew symmetric gradient operator, formally equivalent to `f -> ½(∇f - (∇f)ᵀ)`.
"""
skew_symmetric_gradient(f) = Operation(skew_symmetric_part)(gradient(f))

function evaluate!(cache,::Broadcasting{typeof(skew_symmetric_gradient)},f)
  Broadcasting(Operation(skew_symmetric_part))(Broadcasting(∇)(f))
end

"""
    curl(f)

Abstract curl operator, formally equivalent to
- `f -> ∂₁f₂ - ∂₂f₁` for 2D vector functions, or
- `f -> ∇×f` for 3D vector functions.
"""
curl(f) = Operation(grad2curl)(∇(f))

function evaluate!(cache,::Broadcasting{typeof(curl)},f)
  Broadcasting(Operation(grad2curl))(Broadcasting(∇)(f))
end

"""
    grad2curl(∇f)

Return
- `∇f[1,2] - ∇f[2,1]` for 2×2 input tensor, or
- `VectorValue(∇f[2,3] - ∇f[3,2], ∇f[3,1] - ∇f[1,3], ∇f[1,2] - ∇f[2,1])` for 3×3 input tensor.
"""
function grad2curl(∇u::TensorValue{2})
  ∇u[1,2] - ∇u[2,1]
end

function grad2curl(∇u::TensorValue{3})
  c1 = ∇u[2,3] - ∇u[3,2]
  c2 = ∇u[3,1] - ∇u[1,3]
  c3 = ∇u[1,2] - ∇u[2,1]
  VectorValue(c1,c2,c3)
end

function laplacian end

"""
    const Δ = laplacian

Alias for `laplacian`.
"""
const Δ = laplacian

"""
    laplacian(f)

Abstract laplacian operator, equivalent to `tr(∇∇(f))`.
"""
function laplacian(f)
  # g = gradient(f)
  # divergence(g)
  tr(∇∇(f))
end

"""
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
