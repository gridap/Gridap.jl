# @santiagobadia : Do we need this?
function testitem(af::AbstractArray{<:NewField},x...)
  evaluate(af,x...)
end

@inline function testitem!(cache,af::AbstractArray{<:NewField},x...)
  evaluate!(cache,af,x...)
end


# @inline function gradient(af::AbstractArray{<:NewField})
#   map(∇,af)
# end
