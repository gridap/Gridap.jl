using Gridap.Helpers
using Gridap.FieldValues
using Gridap.Maps
import Gridap: evaluate!, return_size, gradient

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

struct MockGeoMapJaco{D,DD} <: Field{D,TensorValue{D,DD}} end

return_size(this::MockGeoMapJaco, psize::Tuple{Int}) = psize

function evaluate!(
  this::MockGeoMapJaco{D,DD},
  a::AbstractVector{Point{D}},
  v::AbstractVector{TensorValue{D,DD}}) where {D,DD}
  t = TensorValue(2.0,0.0,0.0,2.0)
  for i in eachindex(v)
    v[i] = t
  end
end

function gradient(this::MockGeoMapJaco)
  @notimplemented
end

struct MockGeoMap{D} <: Geomap{D,D} end

function MockGeoMap(D)
  MockGeoMap{D}()
end

return_size(this::MockGeoMap, psize::Tuple{Int}) = psize

function evaluate!(
  this::MockGeoMap{D},
  a::AbstractVector{Point{D}},
  v::AbstractVector{Point{D}}) where D
  v .= 2*a
end

function gradient(this::MockGeoMap{D}) where D
  MockGeoMapJaco{D,D*D}()
end

