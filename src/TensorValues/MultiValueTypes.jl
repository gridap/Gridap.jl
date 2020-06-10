###############################################################
# MultiValue Type
###############################################################

"""
Type representing a multi-dimensional value
"""
abstract type MultiValue{S,T,N,L} <: Number end

@inline Base.Tuple(arg::MultiValue) = arg.data

# Custom type printing

function show(io::IO,v::MultiValue)
  print(io,v.data)
end

function show(io::IO,::MIME"text/plain",v:: MultiValue)
  print(io,typeof(v))
  print(io,v.data)
end

###############################################################
# Other constructors and conversions implemented for more generic types
###############################################################

change_eltype(::Type{<:Number},::Type{T}) where {T} = T
change_eltype(::Number,::Type{T2}) where {T2} = change_eltype(Number,T2)

num_components(::Type{<:Number}) = 1
num_components(::Number) = num_components(Number)

function n_components(a)
  msg = "Function n_components has been removed, use num_components instead"
  error(msg)
end

function data_index(::Type{<:MultiValue},i...)
  @abstractmethod
end
