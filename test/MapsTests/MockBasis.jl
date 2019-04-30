struct MockBasis{D} <: Map{Point{D},1,Point{D},2}
  val::Point{D}
  dim::Int
end

function evaluate!(
  this::MockBasis{D},
  points::AbstractVector{Point{D}},
  v::AbstractArray{Point{D},2}) where {D}
  for j in 1:this.dim
    for (i,pi) in enumerate(points)
      v[j,i] = j*pi+this.val
    end
  end
end

return_size(this::MockBasis, psize::Tuple{Int}) = (this.dim, psize...)

function gradient(this::MockBasis{D}) where D
  MockBasis{D}(zeros(Point{D}), this.dim)
end
