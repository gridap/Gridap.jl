# +++ Important! +++
# Default implementation of broadcasting
#
# In general the broadcasted function can be anything including, gradient, ∘, etc.
# The only requirement is that it is defined for Field objects.
#
# Since overloading the Julia broadcast machinery is challenging,
# we overload this function signature for optimizations, e.g.
#
# function evaluate!(cache,b::BroadcastMapping{typeof(gradient)},g::AbstractArray{<:Field}})
#
# implements the efficient version. From the user perspective all these do conceptually the same:
#
#  ∇.(i_to_f)
#  broadcast(∇,i_to_f)
#  evaluate(BroadcastMapping(∇),i_to_f) or simply BroadcastMapping(∇)(i_to_f)
#
#  but the BroadcastMapping is the one implemented efficiently.
#
function evaluate!(cache,b::BroadcastMapping,args::Union{Field,AbstractArray{<:Field}}...)
  broadcast(b.f,args...)
end
