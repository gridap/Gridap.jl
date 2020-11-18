
"""
    divergence(f)
"""
divergence(f) = Operation(tr)(∇(f))

function evaluate!(cache,::Broadcasting{typeof(divergence)},f)
  Broadcasting(Operation(tr))(Broadcasting(∇)(f))
end

function symmetric_gradient end

"""
    symmetric_gradient(f)
"""
symmetric_gradient(f) = Operation(symmetric_part)(gradient(f))

function evaluate!(cache,::Broadcasting{typeof(symmetric_gradient)},f)
  Broadcasting(Operation(symmetric_part))(Broadcasting(∇)(f))
end

"""
    const ε = symmetric_gradient

Alias for the symmetric gradient
"""
const ε = symmetric_gradient

"""
    curl(f)
"""
curl(f) = Operation(grad2curl)(∇(f))

function evaluate!(cache,::Broadcasting{typeof(curl)},f)
  Broadcasting(Operation(grad2curl))(Broadcasting(∇)(f))
end

"""
    grad2curl(∇f)
"""
@inline function grad2curl(∇u::TensorValue{2})
  ∇u[1,2] - ∇u[2,1]
end

@inline function grad2curl(∇u::TensorValue{3})
  c1 = ∇u[2,3] - ∇u[3,2]
  c2 = ∇u[3,1] - ∇u[1,3]
  c3 = ∇u[1,2] - ∇u[2,1]
  VectorValue(c1,c2,c3)
end

function laplacian end

"""
    const Δ = laplacian

Alias for the `laplacian` function
"""
const Δ = laplacian

"""
    laplacian(f)
"""
function laplacian(f)
  g = gradient(f)
  divergence(g)
end

"""
    ∇⋅f

Equivalent to

    divergence(f)
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

Equivalent to

    gradient(f)
"""
outer(::typeof(∇),f::Field) = gradient(f)
outer(::typeof(∇),f::Function) = gradient(f)

"""
    outer(f,∇)

Equivalent to

    transpose(gradient(f))
"""
outer(f::Field,::typeof(∇)) = transpose(gradient(f))
outer(f::Function,::typeof(∇)) = transpose(gradient(f))

"""
    cross(∇,f)

Equivalent to
    
    curl(f)
"""
cross(::typeof(∇),f::Field) = curl(f)
cross(::typeof(∇),f::Function) = curl(f)

