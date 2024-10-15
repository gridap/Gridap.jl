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

"""
    indep_components_names(::MultiValue)

Returns an array of strings containing the component labels in the order they are stored internally, consistently with _prepare_data(::Multivalue)

If all dimensions of the tensor shape S are smaller than 3, the components should be named with letters "X","Y" and "Z" similarly to the automatic naming of Paraview. Else, if max(S)>3, they are labeled from "1" to "\$dim".
"""
function indep_components_names(::Type{MultiValue{S,T,N,L}}) where {S,T,N,L}
  return ["$i" for i in 1:L]
end

