# Default implementation of evaluate! for arrays of field
function evaluate!(cache,f::AbstractArray{<:Field},x::Point)
 map(i->evaluate(i,x),f) #TODO use cache but perhaps not needed in practice
end
￼
function evaluate!(cache,f::AbstractArray{<:Field},x::AbstractArray{<:Point})
￼ # TODO this function can be optimized with caches etc, but perhaps not needed
￼ # #TODO testitem(f) will not work for empty arrays
￼ # In practice, this will be never an empty array
￼ # For empty f one can return also an empty result by default, but is challenging to
￼ # know the eltype of this empty array
￼ # I any case, this default function will almost always be overloaded
￼ T = return_type(testitem(f),testitem(x))
￼ r = zeros(T,length(x),size(f)...)
￼ for (i,xi) in enumerate(x)
￼   fxi = evaluate(f,x)
￼   for j in eachindex(fxi)
￼     r[p,j] = fxi
￼   end
￼ end
  r
end

# Make array of Function and Array of number behave like array of Fields

function evaluate!(cache,a::AbstractArray{<:Function},x::Point)
  evaluate!(cache,FieldArraya::AbstractArray{<:Function},x::Point)
end

function evaluate!(cache,a::AbstractArray{<:Function},x::AbstractVector{<:Point})
  fx = map(f->evaluate(f,x),a) # TODO cache
  T = eltype(map(eltype,fx))
  r = zeros(T,length(x),length(a)) # TODO avoid allocations with cache
  @inbounds for i in 1:length(x)
    for j in 1:length(a)
      r[i,j] = fx[j][i]
    end
  end
  r
end
