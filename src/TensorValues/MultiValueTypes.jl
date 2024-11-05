###############################################################
# MultiValue Type
###############################################################

"""
    MultiValue{S,T,N,L} <: Number

Abstract type representing a multi-dimensional number value. The parameters are analog to that of StaticArrays.jl:
- `S` is a Tuple type holding the size of the tensor, e.g. Tuple{3} for a 3d vector or Tuple{2,4} for a 2 rows and 4 columns tensor,
- `T` is the type of the scalar components, should be subtype of `Number`,
- `N` is the order of the tensor, the length of `S`,
- `L` is the number of components stored internally.

`MultiValue`s are immutable. See [`TensorValues`](@ref) for more details on usage.
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

"""
    change_eltype(m::Number,::Type{T2})
    change_eltype(M::Type{<:Number},::Type{T2})

For multivalues, returns `M` or `typeof(m)` but with the component type (`MultiValue`'s parametric type `T`) changed to `T2`.

For scalars (or any non MultiValue number), `change_eltype` returns T2.
"""
change_eltype(::Type{<:Number},::Type{T}) where {T} = T
change_eltype(::Number,::Type{T2}) where {T2} = change_eltype(Number,T2)


"""
    Mutable(T::Type{<:MultiValue}) -> ::Type{<:MArray}
    Mutable(a::MultiValue)

Return the concrete `MArray` type (defined by `StaticArrays.jl`) corresponding
to the `MultiValue` type T or array size and type of `a`.

See also [`mutable`](@ref).
"""
Mutable(::Type{MultiValue}) = @abstractmethod
Mutable(::MultiValue) = Mutable(MultiValue)

"""
    mutable(a::MultiValue)

Converts `a` into an array of type `MArray` defined by `StaticArrays.jl`.

See also [`Mutable`](@ref).
"""
mutable(a::MultiValue) = @abstractmethod

"""
    num_components(::Type{<:Number})
    num_components(a::Number)

Total number of components of a `Number` or `MultiValue`, that is 1 for scalars
and the product of the size dimensions for a `MultiValue`. This is the same as `length`.
"""
num_components(::Type{<:Number}) = 1
num_components(::Number) = num_components(Number)
num_components(T::Type{<:MultiValue}) = @unreachable "$T type is too abstract to count its components, provide a (parametric) concrete type"

"""
    num_indep_components(::Type{<:Number})
    num_indep_components(a::Number)

Number of independant components of a `Number`, that is `num_components`
minus the number of components determined from others by symmetries or constraints.

For example, a `TensorValue{3,3}` has 9 independant components, a `SymTensorValue{3}`
has 6 and a `SymTracelessTensorValue{3}` has 5. But they all have 9 (non independant) components.
"""
num_indep_components(::Type{T}) where T<:Number = num_components(T)
num_indep_components(::T) where T<:Number = num_indep_components(T)

function n_components(a)
  msg = "Function n_components has been removed, use num_components instead"
  error(msg)
end

# This should probably not be exported, as (accessing) the data field of
# MultiValue is not a public api
"""
Transforms Cartesian indices to linear indices that index `MultiValue`'s private internal storage, this should'nt be used.
"""
function data_index(::Type{<:MultiValue},i...)
  @abstractmethod
end

# The order of export of components is that of their position in the .data
# field, but the actual method "choosing" the export order is
# Gridap.Visualization._prepare_data(::Multivalue).
"""
    indep_comp_getindex(a::Number,i)

Get the `i`th independent component of `a`. It only differs from `getindex(a,i)`
when the components of `a` are interdependant, see [`num_indep_components`](@ref).
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

"""
    indep_components_names(::MultiValue)

Return an array of strings containing the component labels in the order they
are exported in VTK file.

If all dimensions of the tensor shape S are smaller than 3, the components
are named with letters "X","Y" and "Z" similarly to the automatic naming
of Paraview. Else, if max(S)>3, they are labeled by integers starting from "1".
"""
function indep_components_names(::Type{MultiValue{S,T,N,L}}) where {S,T,N,L}
  return ["$i" for i in 1:L]
end
