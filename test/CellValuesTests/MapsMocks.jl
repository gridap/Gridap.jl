module MapsMocks

using Gridap
using TensorValues

export MockMap
export TestMap

import Gridap: evaluate!
import Gridap: return_size

struct MockMap{P} <: Map{VectorValue,1,P,1}
  val::P
end

function evaluate!(
  this::MockMap{P}, points::AbstractVector{<:VectorValue}, v::AbstractVector{P}) where P
  for (i,qi) in enumerate(points)
    v[i] = qi+this.val
  end
end

return_size(::MockMap, psize::Tuple{Int}) = psize

struct TestMap{P<:VectorValue} <: Map{VectorValue,1,P,2}
  val::P
  dim::Int
end

function evaluate!(
  this::TestMap{P}, points::AbstractVector{<:VectorValue}, v::AbstractMatrix{P}) where P
  for j in 1:this.dim
    for (i,qi) in enumerate(points)
      v[j,i] = j*qi+this.val
    end
  end
end

return_size(this::TestMap, psize::Tuple{Int}) = (this.dim, psize[1])

end # module
