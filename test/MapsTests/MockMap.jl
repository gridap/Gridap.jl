using Numa.Helpers
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

struct TestMap{T,N,D} <: Map{Point{D},1,T,N}
  val::Array{T,N}
end

function TestMap(val::Array{T,N},D) where {T,N}
  TestMap{T,N,D}(val)
end

function evaluate!(
  this::TestMap{T,N,D},
  a::AbstractVector{Point{D}},
  v::AbstractArray{T,N}) where {T,N,D}
  v .= this.val
end

return_size(this::TestMap, psize::Tuple{Int}) = size(this.val)

function gradient(this::TestMap)
  @notimplemented
end

