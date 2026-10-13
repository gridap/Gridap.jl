###############################################################
# Types
###############################################################

"""
    ScalarValue{T} <: MultiValue{Tuple{},T,0,1}

Type representing an order 0 tensor, that is scalar of type `T`.

The component is reached by indexing with no index, or with the linear index
`1`, and the tensor converts to and from `T`, `StaticArrays.Scalar{T}`
(`SArray{Tuple{},T,0,1}`) and `Array{T,0}`:

    a = ScalarValue(2.0)
    2.0 === a[] === a[1] === a[CartesianIndex()] # true
"""
struct ScalarValue{T} <: MultiValue{Tuple{},T,0,1}
  data::T

  function ScalarValue{T}(data::NTuple{1,T}) where T
    T<:Union{MultiValue,AbstractArray} && @unreachable "The component type of a
      ScalarValue must be neither a MultiValue nor an AbstractArray, got $T; use
      `product_tensor` or `outer` to combine tensors"
    new{T}(only(data))
  end
end

function promote_rule(::Type{ScalarValue{Ta}}, ::Type{ScalarValue{Tb}}) where {Ta,Tb}
  ScalarValue{promote_type(Ta,Tb)}
end

###############################################################
# Constructors (ScalarValue)
###############################################################

ScalarValue(data::NTuple{1,T})      where T        = ScalarValue{T}(data)
ScalarValue{T1}(data::NTuple{1,T2}) where {T1,T2}  = ScalarValue{T1}(NTuple{1,T1}(data))

ScalarValue(data::Tuple)                    = @unreachable "A ScalarValue holds a
  single component, got $(length(data))"
ScalarValue{T}(data::Tuple)         where T = @unreachable "A ScalarValue holds a
  single component, got $(length(data))"

ScalarValue(data::Number...)            = ScalarValue(data)
ScalarValue{T}(data::Number...) where T = ScalarValue{T}(data)

ScalarValue(data::AbstractArray{T,0})          where T = ScalarValue{T}((data[],))
ScalarValue{T}(data::AbstractArray{<:Any,0})   where T = ScalarValue{T}((data[],))

###############################################################
# Conversions (ScalarValue)
###############################################################

convert(::Type{<:SArray{Tuple{},T}}, arg::ScalarValue) where T = SArray{Tuple{},T,0,1}(Tuple(arg))
convert(::Type{<:MArray{Tuple{},T}}, arg::ScalarValue) where T = MArray{Tuple{},T,0,1}(Tuple(arg))
convert(::Type{<:ScalarValue{T}}, arg::ScalarValue)    where T = ScalarValue{T}(arg[])
convert(::Type{<:ScalarValue{T}}, arg::ScalarValue{T}) where T = arg

convert(::Type{<:ScalarValue{T}}, arg::_Scalar) where T = ScalarValue{T}(T(arg))
convert(::Type{T}, arg::ScalarValue) where T<:_Scalar = T(arg[])

change_eltype(::Type{<:ScalarValue},::Type{T}) where T = ScalarValue{T}

function change_variance(::Type{<:ScalarValue}, ::Val{Sv}) where Sv
  _check_change_variance(ScalarValue, Sv)
  ScalarValue
end
function change_variance(::Type{<:ScalarValue{T}}, ::Val{Sv}) where {T,Sv}
  _check_change_variance(ScalarValue, Sv)
  ScalarValue{T}
end

MultiValue(a::StaticArray{Tuple{},T,0}) where T = ScalarValue{T}(a[])

###############################################################
# Introspection (ScalarValue)
###############################################################

@inline Base.Tuple(arg::ScalarValue) = (arg.data,)

# The multiplicative neutral is the bare component, like that of a 0-form: `*`
# between two MultiValues is removed in favour of `⋅` and `⊙`, so `a*one(a)==a`
# only holds for a scalar neutral
one(::Type{<:ScalarValue{T}}) where T = one(T)

@propagate_inbounds _get_data(arg::ScalarValue, i) = arg.data

###############################################################
# Indexing (ScalarValue)
###############################################################

@propagate_inbounds getindex(arg::ScalarValue) = arg.data
@propagate_inbounds getindex(arg::ScalarValue, ::CartesianIndex{0}) = arg.data

@propagate_inbounds function getindex(arg::ScalarValue, i::Integer)
  @boundscheck @check isone(i) "A ScalarValue holds a single component, got index $i"
  arg.data
end

###############################################################
# Custom type printing (ScalarValue)
###############################################################

show(io::IO,arg::ScalarValue) = print(io,Tuple(arg))
show(io::IO,::MIME"text/plain",arg::ScalarValue) = (print(io,typeof(arg)); print(io,Tuple(arg)))
