struct MockBasis{D} <: Map{Point{D},1,Point{D},2}
  val::Point{D}
  dim::Int
end

function evaluate!(
  this::MockBasis{D},
  points::AbstractArray{Point{D},1},
  v::AbstractArray{Point{D},2}) where {D}
  for j in 1:this.dim
    for (i,pi) in enumerate(points)
      v[j,i] = j*pi+this.val
    end
  end
end

function evaluate(
  this::MockBasis{D},
  points::AbstractArray{Point{D},1}) where {D}
  v_size = return_size(this, size(points))
  v = Array{Point{D},2}(undef, v_size)
  evaluate!(this,points,v)
  return v
end

return_size(this::MockBasis, p_size::NTuple{N,Int} where N) = (this.dim, p_size...)

function gradient(this::MockBasis{D}) where D
  MockBasis{D}(zeros(Point{D}), this.dim)
end
