
apply_function(f,x...) = apply_mapping(FunctionMapping(f),x...)

# Optimisations

# function apply_mapping(
  # a::MappedArray{T,N,F,<:Fill{<:Type{<:AppliedField}}} where {T,N,F},x::AbstractArray)
  # v, w = _split(a.f)
  # @notimplementedif ! isa(v,Fill)
  # k = v.value
  # apply(k,w...)
# end
