struct MockMap{D} <: Map{Point{D},1,Point{D},1}
  val::Point{D}
end

function evaluate!(
  this::MockMap{D},
  points::AbstractArray{Point{D},1},
  v::AbstractArray{Point{D},1}) where {D}
  l = length(points)
  for (i,pi) in enumerate(points)
    v[i] = pi+this.val
  end
end

function evaluate(
  this::MockMap{D},
  points::AbstractArray{Point{D},1}) where {D}
  v_size = return_size(this, size(points))
  v = Array{Point{D},1}(undef, v_size)
  evaluate!(this,points,v)
  return v
end

return_size(::Map, p_size::NTuple{N,Int} where N) = p_size

function gradient(this::MockMap{D}) where D
  MockMap{D}(zeros(Point{D}))
end
