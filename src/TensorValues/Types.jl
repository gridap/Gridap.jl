
# Types

"""
Type representing a multi-dimensional value
"""
struct MultiValue{S,T,N,L} <: Number
  array::SArray{S,T,N,L}
end

"""
Type representing a second-order tensor
"""
const TensorValue{D,T,L} = MultiValue{Tuple{D,D},T,2,L}

"""
Type representing a first-order tensor
"""
const VectorValue{D,T} = MultiValue{Tuple{D},T,1,D}

# Constructors (MultiValue)

function (::Type{MultiValue{S}})(x::Tuple) where S<:Tuple
  array = SArray{S}(x)
  MultiValue(array)
end

function (::Type{MultiValue{S}})(x::Tuple{}) where S<:Tuple
  s = """
  Unknown element type.

  Provide element type in the corresponding type parameter.
  Examples:
  MultiValue{Tuple{0,0},Int}()
  TensorValue{0,Int}()
  VectorValue{0,Int}()
  """
  error(s)
end

function (::Type{MultiValue{S,T}})(x::Tuple) where {S<:Tuple,T}
  array = SArray{S,T}(x)
  MultiValue(array)
end

function (::Type{MultiValue{S}})(x::Vararg) where S<:Tuple
  MultiValue{S}(x)
end

function (::Type{MultiValue{S,T}})(x::Vararg) where {S<:Tuple,T}
  MultiValue{S,T}(x)
end

function MultiValue(a::StaticArray{S,T}) where {S,T}
  MultiValue{S,T}(a.data)
end

# Constructors (TensorValue)

function (::Type{TensorValue{D}})(x::Tuple) where D
  S = Tuple{D,D}
  MultiValue{S}(x)
end

function (::Type{TensorValue{0}})()
  S = Tuple{0,0}
  MultiValue{S}()
end

function (::Type{TensorValue{D}})(x::Vararg) where D
  TensorValue{D}(x)
end

function (::Type{TensorValue{D,T}})(x::Tuple) where {D,T}
  S = Tuple{D,D}
  MultiValue{S,T}(x)
end

function (::Type{TensorValue{D,T}})(x::Vararg) where {D,T}
  TensorValue{D,T}(x)
end

@generated function TensorValue(arg::NTuple{DD,T}) where {T,DD}
  SQ = sqrt(DD)
  D = ceil(Int,SQ)
  @assert D == SQ
  :( TensorValue{$D,T}(arg)  )
end

function TensorValue(args::Vararg)
  TensorValue(args)
end

function TensorValue()
  S = Tuple{0,0}
  MultiValue{S}()
end

function TensorValue(a::StaticArray)
  TensorValue(a.data)
end

"""
"""
@generated function diagonal_tensor(v::VectorValue{D,T}) where {D,T}
  s = ["zero(T), " for i in 1:(D*D)]
  for i in 1:D
    d = D*(i-1)+i
    s[d] = "v.array[$i],"
  end
  str = join(s)
  Meta.parse("TensorValue(($str))")
end

# Constructors (VectorValue)

function (::Type{VectorValue{D}})(x::Tuple) where D
  S = Tuple{D}
  MultiValue{S}(x)
end

function (::Type{VectorValue{D}})(x::Vararg) where D
  VectorValue{D}(x)
end

function (::Type{VectorValue{D,T}})() where {D,T}
  S = Tuple{D}
  MultiValue{S,T}()
end

function VectorValue(arg::NTuple{D,T}) where {D,T}
  VectorValue{D,T}(arg)
end

function (::Type{VectorValue{D,T}})(x::Vararg{Number,D}) where {T,D}
  VectorValue{D,T}(x)
end

function VectorValue(args::Vararg)
  VectorValue(args)
end

function VectorValue()
  S = Tuple{0}
  MultiValue{S}()
end

function VectorValue(a::StaticArray)
  VectorValue(a.data)
end

function VectorValue(a::SVector)
  MultiValue(a)
end

function VectorValue(a::MVector)
  MultiValue(a)
end



# Initializers

function zero(::Type{<:MultiValue{S,T,N,L}}) where {S,T,N,L}
  z = zero(SArray{S,T,N,L})
  MultiValue{S,T,N,L}(z)
end

function zero(::MultiValue{S,T,N,L}) where {S,T,N,L}
  zero(MultiValue{S,T,N,L})
end

function one(::Type{<:MultiValue{S,T,N,L}}) where {S,T,N,L}
  z = one(SArray{S,T,N,L})
  MultiValue{S,T,N,L}(z)
end

function one(::MultiValue{S,T,N,L}) where {S,T,N,L}
  one(MultiValue{S,T,N,L})
end

# Conversions

function convert(::Type{<:MultiValue{S,T,N,L}},a::StaticArray{S,T,N}) where {S,T,N,L}
  MultiValue(a)
end

function convert(
  ::Type{<:MultiValue{S,T,N,L}},a::AbstractArray{R,N}) where {S,T,N,L,R}
  b = convert(SArray{S,T,N,L},a)
  MultiValue(b)
end

function convert(::Type{<:MultiValue{S,T,N,L}},a::NTuple{L,R}) where {S,T,N,L,R}
  MultiValue(SArray{S,T}(a))
end

function convert(::Type{<:MultiValue{S,T,N,L}},a::MultiValue{S,Ta,N,L} where Ta) where {S,T,N,L}
  b = convert(SArray{S,T,N,L},a.array)
  MultiValue(b)
end

# Misc operations on the type itself

length(::Type{<: MultiValue{S,T,N,L} where {S,T,N}}  ) where L = L

function size(::Type{MultiValue{S,T,N,L}}) where {S,T,N,L}
  A = SArray{S,T,N,L}
  size(A)
end

function size(::Type{<:MultiValue{S}}) where S
  _s(Size(S))
end

@pure _s(s::Size{T}) where T = T

"""
    n_components(::Type)

Returns the number of components stored in the given type.
Implemented for types `<:Real` and `<:MultiValue`.
Also available for instances of these types.
"""
n_components(::Type{<: MultiValue{S,T,N,L} where {S,T,N}}  ) where L = L
n_components(a::T) where T<:MultiValue = n_components(T)

n_components(::Type{<:Real}) = 1
n_components(a::T) where T<:Real = 1


# Custom type printing

function show(io::IO,v::MultiValue)
  print(io,v.array.data)
end

function show(io::IO,::MIME"text/plain",v::MultiValue)
  print(io,typeof(v))
  print(io,v.array.data)
end

# Misc

"""
"""
mutable(::Type{MultiValue{S,T,N,L}}) where {S,T,N,L} = MArray{S,T,N,L}

mutable(::MultiValue{S,T,N,L}) where {S,T,N,L} = MArray{S,T,N,L}

"""
"""
function change_eltype(::Type{MultiValue{S,T,N,L}},::Type{E}) where {S,T,N,L,E}
  MultiValue{S,E,N,L}
end

change_eltype(a::T,::Type{E}) where {T<:MultiValue,E} = change_eltype(T,E)

change_eltype(::Type{<:Real},::Type{E}) where E = E

change_eltype(a::T,::Type{E}) where {T<:Real,E} = change_eltype(T,E)

@inline Tuple(a::MultiValue) = a.array.data

