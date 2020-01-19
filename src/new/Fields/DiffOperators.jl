
"""
    divergence(f)
"""
divergence(f) = tr(gradient(f))

function symmetric_gradient end

"""
    symmetric_gradient(f)
"""
symmetric_gradient(f) = symmetic_part(gradient(f))

"""
    const ε = symmetric_gradient

Alias for the symmetric gradient
"""
const ε = symmetric_gradient

"""
    curl(f)
"""
curl(f) = grad2curl(gradient(f))

"""
    grad2curl(∇f)
"""
function grad2curl(f)
  @abstractmethod
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
  h = gradient(g)
  tr(h)
end

"""
    ∇*f

Equivalent to

    divergence(f)
"""
(*)(::typeof(∇),f) = divergence(f)
(*)(::typeof(∇),f::GridapType) = divergence(f)

"""
    outer(∇,f)

Equivalent to

    gradient(f)
"""
outer(::typeof(∇),f) = gradient(f)
outer(::typeof(∇),f::GridapType) = gradient(f)

"""
    outer(f,∇)

Equivalent to

    transpose(gradient(f))
"""
outer(f,::typeof(∇)) = transpose(gradient(f))
outer(f::GridapType,::typeof(∇)) = transpose(gradient(f))

"""
    cross(∇,f)

Equivalent to
    
    curl(f)
"""
cross(::typeof(∇),f) = curl(f)
cross(::typeof(∇),f::GridapType) = curl(f)

# Helpers

grad2curl(f::Field) = apply_kernel_to_field(bcast(_curl_kernel),f)

grad2curl(f::AbstractArray{<:Field}) = apply_to_field_array(bcast(_curl_kernel),f)

function _curl_kernel(∇u::TensorValue{2})
  ∇u[1,2] - ∇u[2,1]
end

function _curl_kernel(∇u::TensorValue{3})
  c1 = ∇u[2,3] - ∇u[3,2]
  c2 = ∇u[3,1] - ∇u[1,3]
  c3 = ∇u[1,2] - ∇u[2,1]
  VectorValue(c1,c2,c3)
end

