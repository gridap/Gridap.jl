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
num_components(T::Type{<:MultiValue}) = @unreachable "$T type is too abstract to count its components, provide a (parametric) concrete type"

"""
Number of independant components, that is `num_components(::Type{T})` minus the
number of components determined from others by symmetries or constraints.
"""
num_indep_components(::Type{T}) where T<:Number = num_components(T)
num_indep_components(::T) where T<:Number = num_indep_components(T)

function n_components(a)
  msg = "Function n_components has been removed, use num_components instead"
  error(msg)
end

# This should probably not be exported, as (accessing) the data field of
# MultiValue is not a public api
function data_index(::Type{<:MultiValue},i...)
  @abstractmethod
end

"""
    indep_comp_getindex(a::Number,i)

Get the ith independent component of `a`. It only differs from `getindex(a,i)`
when the components of `a` are linked, see [`num_indep_components`](@ref), and
`i` should be in `1:num_indep_components(a)`.
"""
function indep_comp_getindex(a::Number,i)
  @check 1 <= i <= num_indep_components(Number)
  a[i]
end

function indep_comp_getindex(a::T,i) where {T<:MultiValue}
  @check 1 <= i <= num_indep_components(T)
  _get_data(a,i)
end

# abstraction of Multivalue data access in case subtypes of MultiValue don't
# store its data in a data field
function _get_data(a::MultiValue,i)
  a.data[i]
end
