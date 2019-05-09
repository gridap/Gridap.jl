using Numa.FieldValues
using Numa.Maps
import Numa: evaluate!, return_size, gradient

struct MockMap{D} <: Map{Point{D},1,Point{D},1}
  val::Point{D}
end

function evaluate!(
  this::MockMap{D},
  points::AbstractVector{Point{D}},
  v::AbstractVector{Point{D}}) where {D}
  for (i,qi) in enumerate(points)
    v[i] = qi+this.val
  end
end

return_size(::MockMap, psize::Tuple{Int}) = psize

function gradient(this::MockMap{D}) where D
  MockMap{D}(zeros(Point{D}))
end
