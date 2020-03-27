###############################################################
# Types
###############################################################

"""
Type representing a multi-dimensional value
"""
abstract type MultiValue{S,T,N,L} <: Number end

"""
Type representing a first-order tensor
"""
struct VectorValue{D,T} <: MultiValue{Tuple{D},T,1,D}
    data::NTuple{D,T}
    function VectorValue{D,T}(data::NTuple{D,T}) where {D,T}
        new{D,T}(data)
    end
end

"""
Type representing a second-order tensor
"""
struct TensorValue{D1,D2,T,L} <: MultiValue{Tuple{D1,D2},T,2,L}
    data::NTuple{L,T}
    function TensorValue{D1,D2,T}(data::NTuple{L,T}) where {D1,D2,T,L}
        @assert L == D1*D2
        new{D1,D2,T,L}(data)
    end
end

###############################################################
# Constructors (VectorValue)
###############################################################

# VectorValue no-arguments constructor

function VectorValue() 
    VectorValue{0}()
end

function VectorValue{0}() 
    VectorValue{0,Int}()
end

function VectorValue{0,T}() where {T}
    VectorValue{0,T}(NTuple{0,T}())
end

# VectorValue empty tuple constructor

function VectorValue(data::NTuple{0})
    VectorValue{0,Int}(data)
end

function VectorValue{0}(data::NTuple{0})
    VectorValue{0,Int}(data)
end

# VectorValue single NTuple argument constructor

function VectorValue(data::NTuple{D,T}) where {D,T}
    VectorValue{D,T}(data)
end

function VectorValue{D}(data::NTuple{D,T}) where {D,T}
    VectorValue{D,T}(data)
end

function VectorValue{D,T1}(data::NTuple{D,T2}) where {D,T1,T2}
    VectorValue{D,T1}(NTuple{D,T1}(data))
end

# VectorValue Vararg constructor

function VectorValue(data::Vararg)
    VectorValue(NTuple{length(data)}(data))
end

function VectorValue{D}(data::Vararg) where {D}
    VectorValue{D}(NTuple{D}(data))
end

function VectorValue{D,T}(data::Vararg) where {D,T}
    VectorValue{D,T}(NTuple{D,T}(data))
end

# VectorValue single MVector argument constructor

function VectorValue(data::MVector{D,T}) where {D,T}
    VectorValue{D,T}(NTuple{D,T}(data))
end

function VectorValue{D}(data::MVector{D,T}) where {D,T}
    VectorValue{D,T}(NTuple{D,T}(data))
end

function VectorValue{D,T1}(data::MVector{D,T2}) where {D,T1,T2}
    VectorValue{D,T1}(NTuple{D,T1}(data))
end

# VectorValue single MVector argument constructor

function VectorValue(data::SVector{D,T}) where {D,T}
    VectorValue{D,T}(NTuple{D,T}(data))
end

function VectorValue{D}(data::SVector{D,T}) where {D,T}
    VectorValue{D,T}(NTuple{D,T}(data))
end

function VectorValue{D,T1}(data::SVector{D,T2}) where {D,T1,T2}
    VectorValue{D,T1}(NTuple{D,T1}(data))
end

# VectorValue single AbstractArray argument constructor

function VectorValue(data::AbstractArray{T}) where {T}
    D = length(data)
    VectorValue{D,T}(NTuple{D,T}(data))
end

function VectorValue{D}(data::AbstractArray{T}) where {D,T}
    VectorValue{D,T}(NTuple{D,T}(data))
end

function VectorValue{D,T1}(data::AbstractArray{T2}) where {D,T1,T2}
    VectorValue{D,T1}(NTuple{D,T1}(data))
end

###############################################################
# Constructors (TensorValue)
###############################################################

# TensorValue no-arguments constructor

function TensorValue() 
    TensorValue{0,0}()
end

function TensorValue{0,0}() 
    TensorValue{0,0,Int}(NTuple{0,Int}())
end

function TensorValue{0,0,T}() where {T}
    TensorValue{0,0,T}(NTuple{0,T}())
end

# TensorValue empty tuple constructor

function TensorValue(data::NTuple{0})
    TensorValue{0,0,Int}(data)
end

function TensorValue{0,0}(data::NTuple{0})
    TensorValue{0,0,Int}(data)
end

# TensorValue single NTuple argument constructor

function TensorValue(data::NTuple{L,T}) where {L,T}
    D=Int(sqrt(L))
    TensorValue{D,D,T}(data)
end

function TensorValue{D}(data::NTuple{L,T}) where {D,L,T}
    TensorValue{D,D,T}(data)
end

function TensorValue{D1,D2}(data::NTuple{L,T}) where {D1,D2,L,T}
    TensorValue{D1,D2,T}(data)
end

function TensorValue{D1,D2,T1}(data::NTuple{L,T2}) where {D1,D2,L,T1,T2}
    TensorValue{D1,D2,T1}(NTuple{L,T1}(data))
end

function TensorValue{D1,D2,T1,L}(data::NTuple{L,T2}) where {D1,D2,L,T1,T2}
    TensorValue{D1,D2,T1}(NTuple{L,T1}(data))
end

# TensorValue Vararg constructor

function TensorValue(data::Vararg)
    TensorValue(NTuple{length(data)}(data))
end

function TensorValue{D}(data::Vararg) where {D}
    TensorValue{D,D}(NTuple{D*D}(data))
end

function TensorValue{D1,D2}(data::Vararg) where {D1,D2}
    TensorValue{D1,D2}(NTuple{D1*D2}(data))
end

function TensorValue{D1,D2,T}(data::Vararg) where {D1,D2,T}
    TensorValue{D1,D2,T}(NTuple{D1*D2,T}(data))
end

# TensorValue single MVector argument constructor

function TensorValue(data::MVector{L,T}) where {L,T}
    D=Int(sqrt(L))
    TensorValue{D,D,T}(NTuple{L,T}(data))
end

function TensorValue{D1,D2}(data::MVector{L,T}) where {D1,D2,T,L}
    TensorValue{D1,D2,T}(NTuple{L,T}(data))
end

function TensorValue{D1,D2,T1}(data::MVector{L,T2}) where {D1,D2,T1,T2,L}
    TensorValue{D1,D2,T1}(NTuple{L,T1}(data))
end

# TensorValue single SVector argument constructor

function TensorValue(data::SVector{L,T}) where {L,T}
    D=Int(sqrt(L))
    TensorValue{D,D,T}(NTuple{L,T}(data))
end

function TensorValue{D1,D2}(data::SVector{L,T}) where {D1,D2,T,L}
    TensorValue{D1,D2,T}(NTuple{L,T}(data))
end

function TensorValue{D1,D2,T1}(data::SVector{L,T2}) where {D1,D2,T1,T2,L}
    TensorValue{D1,D2,T1}(NTuple{L,T1}(data))
end

# TensorValue single MMAtrix argument constructor

function TensorValue(data::StaticArrays.MMatrix{D1,D2,T,L}) where {D1,D2,T,L}
    D=Int(sqrt(L))
    TensorValue{D,D,T}(NTuple{L,T}(data))
end

function TensorValue{D1,D2}(data::StaticArrays.MMatrix{D1,D2,T,L}) where {D1,D2,T,L}
    TensorValue{D1,D2,T}(NTuple{L,T}(data))
end

function TensorValue{D1,D2,T1}(data::StaticArrays.MMatrix{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    TensorValue{D1,D2,T1}(NTuple{L,T1}(data))
end

function TensorValue{D1,D2,T1,L}(data::StaticArrays.MMatrix{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    TensorValue{D1,D2,T1}(NTuple{L,T1}(data))
end

# TensorValue single SMatrix argument constructor

function TensorValue(data::StaticArrays.SMatrix{D1,D2,T,L}) where {D1,D2,T,L}
    D=Int(sqrt(L))
    TensorValue{D,D,T}(NTuple{L,T}(data))
end

function TensorValue{D1,D2}(data::StaticArrays.SMatrix{D1,D2,T,L}) where {D1,D2,T,L}
    TensorValue{D1,D2,T}(NTuple{L,T}(data))
end

function TensorValue{D1,D2,T1}(data::StaticArrays.SMatrix{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    TensorValue{D1,D2,T1}(NTuple{L,T1}(data))
end

function TensorValue{D1,D2,T1,L}(data::StaticArrays.SMatrix{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    TensorValue{D1,D2,T1}(NTuple{L,T1}(data))
end

# TensorValue single AbstractArray argument constructor

function TensorValue(data::AbstractArray{T}) where {T}
    L = length(data)
    TensorValue(NTuple{L,T}(data))
end

function TensorValue{D1,D2}(data::AbstractArray{T}) where {D1,D2,T}
    TensorValue{D1,D2,T}(NTuple{D1*D2,T}(data))
end

function TensorValue{D1,D2,T1}(data::AbstractArray{T2}) where {D1,D2,T1,T2}
    TensorValue{D1,D2,T1}(NTuple{D1*D2,T1}(data))
end

###############################################################
# Conversions (VectorValue)
###############################################################

function convert(::T, arg::NTuple{D}) where {D,T<:VectorValue}
@show T
    T(arg)
end

function convert(::Type{<:NTuple}, arg::VectorValue{D,T}) where {D,T}
    NTuple{D,T}(arg.data)
end

function convert(::Type{<:NTuple{D}}, arg::VectorValue{D,T}) where {D,T}
    NTuple{D,T}(arg.data)
end

function convert(::Type{<:NTuple{D,T1}}, arg::VectorValue{D,T2}) where {D,T1,T2}
    NTuple{D,T1}(arg.data)
end

function convert(::Type{<:VectorValue}, arg::SVector{D,T}) where {D,T}
    VectorValue{D,T}(arg)
end

function convert(::Type{<:VectorValue{D}}, arg::SVector{D,T}) where {D,T}
    VectorValue{D,T}(arg)
end

function convert(::Type{<:VectorValue{D,T1}}, arg::SVector{D,T2}) where {D,T1,T2}
    VectorValue{D,T1}(arg)
end

function convert(::Type{<:VectorValue}, arg::MVector{D,T}) where {D,T}
    VectorValue{D,T}(arg)
end

function convert(::Type{<:VectorValue{D}}, arg::MVector{D,T}) where {D,T}
    VectorValue{D,T}(arg)
end

function convert(::Type{<:VectorValue{D,T1}}, arg::MVector{D,T2}) where {D,T1,T2}
    VectorValue{D,T1}(arg)
end

function convert(::Type{<:SVector}, arg::VectorValue{D,T}) where {D,T}
    SVector{D,T}(arg.data)
end

function convert(::Type{<:SVector{D}}, arg::VectorValue{D,T}) where {D,T}
    SVector{D,T}(arg.data)
end

function convert(::Type{<:SVector{D,T1}}, arg::VectorValue{D,T2}) where {D,T1,T2}
    SVector{D,T1}(arg.data)
end

function convert(::Type{<:MVector}, arg::VectorValue{D,T}) where {D,T}
    MVector{D,T}(arg.data)
end

function convert(::Type{<:MVector{D}}, arg::VectorValue{D,T}) where {D,T}
    MVector{D,T}(arg.data)
end

function convert(::Type{<:MVector{D,T1}}, arg::VectorValue{D,T2}) where {D,T1,T2}
    MVector{D,T1}(arg.data)
end

function convert(AAT::Type{<:AbstractArray{T1,N}}, arg::VectorValue{D,T2}) where {D,T1,T2,N}
    AAT{T1,N}(collect(T1,arg.data))
end

function convert(::Type{<:VectorValue}, arg::AbstractArray{T,1}) where {T}
    D =length(arg)
    VectorValue{D,T}(NTuple{D,T}(arg))
end

function convert(::Type{<:VectorValue{D}}, arg::AbstractArray{T,1}) where {D,T}
    VectorValue{D,T}(NTuple{D,T}(arg))
end

function convert(::Type{<:VectorValue{D,T1}}, arg::AbstractArray{T2,1}) where {D,T1,T2}
    VectorValue{D,T1}(NTuple{D,T1}(arg))
end

function convert(::Type{<:VectorValue}, arg::VectorValue{D,T2}) where {D,T2}
    arg
end

function convert(::Type{<:VectorValue{D}}, arg::VectorValue{D,T2}) where {D,T2}
end

function convert(::Type{<:VectorValue{D,T1}}, arg::VectorValue{D,T2}) where {D,T1,T2}
    T1 == T2 ? arg : convert(VectorValue{D,T1}, arg.data)
end

###############################################################
# Conversions (TensorValue)
###############################################################

function convert(::Type{<:TensorValue}, arg::NTuple{L}) where {L}
    TensorValue(arg)
end

function convert(::Type{<:TensorValue{D1,D2}}, arg::NTuple{L}) where {D1,D2,L}
    TensorValue{D1,D2}(arg)
end

function convert(::Type{<:TensorValue{D1,D2,T1}}, arg::NTuple{L}) where {D1,D2,T1,L}
    TensorValue{D1,D2,T1}(arg)
end

function convert(::Type{<:NTuple}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T2,L}
    NTuple(arg.data)
end

function convert(::Type{<:NTuple{L,T1}}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    NTuple{L,T1}(arg.data)
end

function convert(::Type{<:TensorValue}, arg::StaticArrays.SMatrix{D1,D2,T2,L}) where {D1,D2,T2,L}
    TensorValue{D1,D2,T2}(arg)
end

function convert(::Type{<:TensorValue{D1,D2,T1}}, arg::StaticArrays.SMatrix{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    TensorValue{D1,D2,T1}(arg)
end

function convert(::Type{<:TensorValue}, arg::StaticArrays.MMatrix{D1,D2,T2,L}) where {D1,D2,T2,L}
    TensorValue{D1,D2,T2}(arg)
end

function convert(::Type{<:TensorValue{D1,D2,T1}}, arg::StaticArrays.MMatrix{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    TensorValue{D1,D2,T1}(arg)
end

function convert(::Type{<:StaticArrays.SMatrix}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T2,L}
    StaticArrays.SMatrix{D1,D2,T2,L}(arg.data)
end

function convert(::Type{<:StaticArrays.SMatrix{D1,D2,T1}}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    StaticArrays.SMatrix{D1,D2,T1,L}(arg.data)
end

function convert(::Type{<:StaticArrays.MMatrix}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T2,L}
    StaticArrays.MMatrix{D1,D2,T2,L}(arg.data)
end

function convert(::Type{<:StaticArrays.MMatrix{D1,D2,T1}}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    StaticArrays.MMatrix{D1,D2,T1,L}(arg.data)
end

function convert(::Type{<:TensorValue}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T2,L}
    arg
end

function convert(::Type{<:TensorValue{D1,D2}}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T2,L}
    arg
end

function convert(::Type{<:TensorValue{D1,D2,T1}}, arg::TensorValue{D1,D2,T2,L}) where {D1,D2,T1,T2,L}
    T1 == T2 ? arg : convert(TensorValue{D1,D2,T1,L}, arg.data)
end

###############################################################
# Other constructors (VectorValue)
###############################################################

function zero(::Type{<:VectorValue{D,T}}) where {D,T}
   VectorValue{D,T}(NTuple{D,T}(zeros(T,D)))
end

function one(::Type{<:VectorValue{D,T}}) where {D,T}
    VectorValue{D,T}(NTuple{D,T}(ones(T,D)))
end

function mutable(::Type{VectorValue{D,T}}) where {D,T}
    MVector{D,T}
end

function change_eltype(::Type{VectorValue{D}},::Type{T}) where {D,T}
    VectorValue{D,T}
end

function change_eltype(::Type{VectorValue{D,T1}},::Type{T2}) where {D,T1,T2}
    VectorValue{D,T2}
end

function zero(::IT where {IT<:VectorValue{D,T}}) where {D,T}
   zero(VectorValue{D,T})
end

function one(::IT where {IT<:VectorValue{D,T}}) where {D,T}
    one(VectorValue{D,T})
end

function mutable(::IT where {IT<:VectorValue{D,T}}) where {D,T}
    mutable(VectorValue{D,T})
end

function change_eltype(::IT where {IT<:VectorValue{D,T1}},::Type{T2}) where {D,T1,T2}
    change_eltype(VectorValue{D,T1},T2)
end

function SVector(arg::IT where {IT<:VectorValue{D,T}}) where {D,T}
    SVector{D,T}(arg.data)
end

function SArray(arg::IT where {IT<:VectorValue{D,T}}) where {D,T}
    SVector(arg)
end

###############################################################
# Other constructors (TensorValue)
###############################################################

function zero(::Type{<:TensorValue{D1,D2,T}}) where {D1,D2,T}
    D=D1*D2
    TensorValue{D1,D2,T}(NTuple{D,T}(zeros(T,D)))
end

@generated function one(::Type{<:TensorValue{D1,D2,T}}) where {D1,D2,T}
  str = join(["$i==$j ? one(T) : zero(T), " for i in 1:D1 for j in 1:D2])
  Meta.parse("TensorValue{D1,D2,T}(($str))")
end

function mutable(::Type{<:TensorValue{D1,D2,T}}) where {D1,D2,T}
    MMatrix{D1,D2,T}
end

function change_eltype(::Type{TensorValue{D1,D2,T1,L}},::Type{T2}) where {D1,D2,T1,T2,L}
    TensorValue{D1,D2,T2,L}
end

function zero(::IT where {IT<:TensorValue{D1,D2,T}}) where {D1,D2,T}
    zero(TensorValue{D1,D2,T})
end

function one(::IT where {IT<:TensorValue{D1,D2,T}}) where {D1,D2,T}
    one(TensorValue{D1,D2,T})
end

function mutable(::IT where {IT<:TensorValue{D1,D2,T}}) where {D1,D2,T}
    mutable(TensorValue{D1,D2,T})
end

function change_eltype(::IT where {IT<:TensorValue{D1,D2,T1,L}},::Type{T2}) where {D1,D2,T1,T2,L}
    change_eltype(TensorValue{D1,D2,T1,L},T2)
end

function StaticArrays.SMatrix(arg::IT where {IT<:TensorValue{D1,D2,T,L}}) where {D1,D2,T,L}
    StaticArrays.SMatrix{D1,D2,T,L}(arg.data)
end

function SArray(arg::IT where {IT<:TensorValue{D1,D2,T,L}}) where {D1,D2,T,L}
    StaticArrays.SMatrix(arg)
end




function get_array(arg::T where {T<:VectorValue}) 
    convert(SVector,arg)
end

function get_array(arg::T where {T<:TensorValue})
    convert(SMatrix,arg)
end

function change_eltype(::Type{<:Number},::Type{T}) where {T}
    T
end

function change_eltype(::Number,::Type{T2}) where {T2}
    change_eltype(Number,T2)
end

function eltype(::Type{<:VectorValue{D,T}}) where {D,T} 
    T
end

function eltype(::Type{<:TensorValue{D1,D2,T}}) where {D1,D2,T} 
    T
end

function eltype(arg::IT where {IT<:VectorValue{D,T}}) where {D,T} 
    eltype(VectorValue{D,T})
end

function eltype(arg::IT where {IT<:TensorValue{D1,D2,T}}) where {D1,D2,T} 
    eltype(TensorValue{D1,D2,T})
end

function size(::Type{VectorValue{D}}) where {D} 
    (D,)
end

function size(::Type{VectorValue{D,T}}) where {D,T} 
    (D,)
end

function size(::Type{TensorValue{D}}) where {D} 
    (D,D)
end

function size(::Type{TensorValue{D1,D2}}) where {D1,D2} 
    (D1,D2)
end

function size(::Type{TensorValue{D1,D2,T}}) where {D1,D2,T} 
    (D1,D2)
end

function size(::Type{TensorValue{D1,D2,T,L}}) where {D1,D2,T,L} 
    (D1,D2)
end

function size(arg::IT where {IT<:VectorValue{D,T}}) where {D,T} 
    size(VectorValue{D,T})
end

function size(arg::IT where {IT<:TensorValue{D1,D2,T}}) where {D1,D2,T} 
    size(TensorValue{D1,D2,T})
end

function length(::Type{VectorValue{D}}) where {D} 
    D
end

function length(::Type{VectorValue{D,T}}) where {D,T} 
    D
end

function length(::Type{TensorValue{D}}) where {D} 
    length(TensorValue{D,D})
end

function length(::Type{TensorValue{D1,D2}}) where {D1,D2} 
    D1*D1
end

function length(::Type{TensorValue{D1,D2,T,L}}) where {D1,D2,T,L} 
    L
end

function length(arg::IT where {IT<:VectorValue{D,T}}) where {D,T} 
    length(VectorValue{D,T})
end

function length(arg::IT where {IT<:TensorValue{D1,D2,T,L}}) where {D1,D2,T,L} 
    length(TensorValue{D1,D2,T,L})
end

function n_components(::Type{<:Number}) 
    1
end

function n_components(::Type{VectorValue{D}}) where {D}
    length(VectorValue{D})
end

function n_components(::Type{VectorValue{D,T}}) where {D,T}
    length(VectorValue{D,T})
end

function n_components(::Type{TensorValue{D}}) where {D}
    length(TensorValue{D,D})
end

function n_components(::Type{TensorValue{D1,D2}}) where {D1,D2}
    length(TensorValue{D1,D2})
end

function n_components(::Type{TensorValue{D1,D2,T,L}}) where {D1,D2,T,L}
    length(TensorValue{D1,D2,T,L})
end

function n_components(arg::IT where {IT<:Number})
    n_components(Number)
end


function n_components(arg::IT where {IT<:VectorValue{D,T}}) where {D,T}
    n_components(VectorValue{D,T})
end


function n_components(arg::IT where {IT<:TensorValue{D1,D2,T,L}}) where {D1,D2,T,L}
    n_components(TensorValue{D1,D2,T,L})
end




"""
"""
@generated function diagonal_tensor(v::VectorValue{D,T}) where {D,T}
  s = ["zero(T), " for i in 1:(D*D)]
  for i in 1:D
    d = D*(i-1)+i
    s[d] = "v.data[$i],"
  end
  str = join(s)
  Meta.parse("TensorValue(($str))")
end

# Misc operations on the type itself

@pure _s(s::Size{T}) where T = T

# Custom type printing

function show(io::IO,v::IT where {IT<:MultiValue})
  print(io,v.data)
end

function show(io::IO,::MIME"text/plain",v:: IT where {IT<:MultiValue})
  print(io,typeof(v))
  print(io,v.data)
end

@inline Tuple(arg::IT where {IT<:MultiValue}) = arg.data

