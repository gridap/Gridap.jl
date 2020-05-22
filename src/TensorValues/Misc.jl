###############################################################
# Other constructors and conversions implemented for more generic types
###############################################################

change_eltype(::Type{<:Number},::Type{T}) where {T} = T
change_eltype(::Number,::Type{T2}) where {T2} = change_eltype(Number,T2)

n_components(::Type{<:Number}) = 1
n_components(::Number) = n_components(Number)

###############################################################
# Misc
###############################################################

# Misc operations on the type itself

@pure _s(s::Size{T}) where T = T

