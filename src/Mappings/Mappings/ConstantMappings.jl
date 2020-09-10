# struct ConstantMapping{T} <: Mapping
#   v::T
#   function ConstantMapping{T}(v::Union{Number,AbstractArray{<:Number}}) where {T}
#     new{typeof(v)}(v)
#   end
# end

# ConstantMapping(v::Union{Number,AbstractArray{<:Number}}) = ConstantMapping{typeof(v)}(v)

# constant_mapping(v::Union{Number,AbstractArray{<:Number}}) = ConstantMapping(v)

# Number

# function return_cache(f::ConstantMapping,x...)
#   nx = length(x)
#   c = zeros(typeof(f.v),nx)
#   CachedArray(c)
# end

# function return_cache(f::ConstantMapping{<:AbstractArray},x::AbstractArray)
#   nx = length(x)
#   sv = size(f.v)
#   s = (nx,sv...)
#   c = zeros(eltype(f.v),s)
#   cis = CartesianIndices(f.v)
#   for i in eachindex(x)
#     for ci in cis
#       @inbounds c[i,ci] = f.v[ci]
#     end
#   end
#   c
# end

# function evaluate!(c,f::ConstantMapping,x...)
#   return f.v
# end

# function evaluate!(c,f::ConstantMapping,x::AbstractArray)
  # return c
# end

# gradient(f::ConstantMapping) = GradientConstantMapping(f)

# # @santiagobadia: Think about recursive struct for higher derivatives
# struct GradientConstantMapping <: Mapping
#   m::ConstantMapping
# end

# function return_cache(f::GradientConstantMapping,x::AbstractArray)
#   nx = length(x)
#   if nx > 1
#     sx = size(x[1]) # must be identical for all entries
#     sv = size(f.m.v)
#     s = (nx,sx...,sv...)
#     c = zeros(eltype(f.m.v),s)
#   else
#     c = nothing
#   end
#   return c
# end

# function evaluate!(c,f::GradientConstantMapping,x...)
#   return c
# end

# function evaluate!(c,f::GradientConstantMapping,x::AbstractArray)
#   return c
# end
