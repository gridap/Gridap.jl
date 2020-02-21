
"""
"""
struct MultiFieldArray{T,N,A<:AbstractArray{T,N}} <: GridapType
  blocks::Vector{A}
  coordinates::Vector{NTuple{N,Int}}
  ptrs::Array{Int,N}
  block_size::NTuple{N,Int}

  function MultiFieldArray(
    blocks::Vector{A},
    coordinates::Vector{NTuple{N,Int}}) where {T,N,A<:AbstractArray{T,N}}

    @assert length(blocks) == length(coordinates)
    msg = "Trying to build a MultiFieldArray with repeated blocks"
    @assert _no_has_repeaded_blocks(coordinates) msg
    ptrs = _prepare_ptrs(coordinates)
    block_size = _get_block_size(coordinates)
    new{T,N,A}(blocks,coordinates,ptrs,block_size)
  end
end

function _no_has_repeaded_blocks(coordinates::Vector{NTuple{N,Int}}) where N
  maxblocks = _get_block_size(coordinates)
  touched = zeros(Int,maxblocks)
  for c in coordinates
    touched[c...] += 1
  end
  all( touched .<= 1 )
end

function _prepare_ptrs(coordinates)
  s = _get_block_size(coordinates)
  ptrs = zeros(Int,s)
  for (i,c) in enumerate(coordinates)
    ptrs[c...] = i
  end
  ptrs
end

function _get_block_size(coordinates::Vector{NTuple{N,Int}}) where N
  m = zeros(Int,N)
  for c in coordinates
    for n in 1:N
      m[n] = max(m[n],c[n])
    end
  end
  NTuple{N,Int}(m)
end

"""
"""
function get_block_size(a::MultiFieldArray)
  a.block_size
end

"""
"""
function num_blocks(a::MultiFieldArray)
  s = get_block_size(a)
  prod(s)
end

"""
"""
function num_stored_blocks(a::MultiFieldArray)
  length(a.coordinates)
end

"""
"""
function has_all_blocks(a::MultiFieldArray{T,N}) where {T,N}
  num_blocks(a) == num_stored_blocks(a)
end

function Base.:*(a::MultiFieldArray,b::Number)
  blocks = [ block*b for block in a.blocks ]
  coordinates = a.coordinates
  MultiFieldArray(blocks,coordinates)
end

function Base.:*(b::Number,a::MultiFieldArray)
  blocks = [ b*block for block in a.blocks ]
  coordinates = a.coordinates
  MultiFieldArray(blocks,coordinates)
end

function Base.show(io::IO,a::MultiFieldArray)
  print(io,"MultiFieldArray($(a.blocks),$(a.coordinates))")
end

function Base.show(io::IO,::MIME"text/plain",a::MultiFieldArray)
  println(io,"MultiFieldArray object:")
  cis = CartesianIndices(a.ptrs)
  for ci in cis
    p = a.ptrs[ci]
    if p == 0
      println(io,"Block $(Tuple(ci)) -> Empty")
    else
      println(io,"Block $(Tuple(ci)) -> $(a.blocks[p])")
    end
  end
end

function add_to_array!(a::MultiFieldArray{Ta,N},b::MultiFieldArray{Tb,N},combine=+) where {Ta,Tb,N}
  for k in 1:length(a.blocks)
    ak = a.blocks[k]
    bk = b.blocks[k]
    add_to_array!(ak,bk,combine)
  end
end

function add_to_array!(a::MultiFieldArray,b::Number,combine=+)
  for k in 1:length(a.blocks)
    ak = a.blocks[k]
    add_to_array!(ak,b,combine)
  end
end

function Base.:*(a::MultiFieldArray{Ta,2},b::MultiFieldArray{Tb,1}) where {Ta,Tb}
  @assert num_stored_blocks(a) != 0
  @notimplementedif ! has_all_blocks(b)

  function fun(i,a,b)
    ai = a.blocks[i]
    ci, cj = a.coordinates[i]
    p = b.ptrs[cj]
    bi = b.blocks[p]
    ai*bi
  end

  blocks = [ fun(i,a,b) for i in 1:length(a.blocks)]
  coordinates = [ (c[1],)  for c in a.coordinates]

  data = _merge_repeated_blocks(blocks,coordinates)
  MultiFieldArray(data...)
end

function Base.fill!(a::MultiFieldArray,b)
  for k in 1:length(a.blocks)
    ak = a.blocks[k]
    fill!(ak,b)
  end
  a
end

Base.eltype(::Type{<:MultiFieldArray{T}}) where T = T

Base.eltype(a::MultiFieldArray{T}) where T = T

function _merge_repeated_blocks(blocks,coordinates::Vector{NTuple{N,Int}}) where N
  @assert length(blocks) == length(coordinates)
  s = _get_block_size(coordinates)
  ptrs = zeros(Int,s)
  A = eltype(blocks)
  _blocks = A[]
  _coords = NTuple{N,Int}[]
  q = 1
  for b in 1:length(blocks)
    c = coordinates[b]
    block = blocks[b]
    p = ptrs[c...]
    if p == 0
      push!(_blocks,block)
      push!(_coords,c)
      ptrs[c...] = q
      q += 1
    else
      add_to_array!(_blocks[p],block)
    end
  end
  (_blocks,_coords)
end

function mul!(c::MultiFieldArray{Tc,1},a::MultiFieldArray{Ta,2},b::MultiFieldArray{Tb,1}) where {Tc,Ta,Tb}
  for ci in c.blocks
    fill!(ci,zero(Tc))
  end
  for k in 1:length(a.blocks)
    ak = a.blocks[k]
    ci, cj = a.coordinates[k]
    p = b.ptrs[cj]
    bk = b.blocks[p]
    q = c.ptrs[ci]
    ck = c.blocks[q]
    matvec_muladd!(ck,ak,bk)
  end
end

function CachedMultiFieldArray(a::MultiFieldArray)
  blocks = [ CachedArray(b) for b in a.blocks ]
  coordinates = a.coordinates
  MultiFieldArray(blocks,coordinates)
end

function _resize_for_mul!(
  c::MultiFieldArray{Tc,1},a::MultiFieldArray{Ta,2},b::MultiFieldArray{Tb,1}) where {Tc,Ta,Tb}
  for k in 1:length(a.blocks)
    ak = a.blocks[k]
    ci, cj = a.coordinates[k]
    q = c.ptrs[ci]
    ck = c.blocks[q]
    setsize!(ck,(size(ak,1),))
  end
end

function _move_cached_arrays!(r::MultiFieldArray,c::MultiFieldArray)
  for  k in 1:length(c.blocks)
    ck = c.blocks[k]
    r.blocks[k] = ck.array
  end
end

function Base.getindex(a::MultiFieldArray{T,N},I::Vararg{Int,N}) where {T,N}
  p = a.ptrs[I...]
  @assert p > 0 "You are attempting to access a block that is not stored"
  a.blocks[p]
end

function Base.getindex(a::MultiFieldArray,i::Integer)
  p = a.ptrs[i]
  @assert p > 0 "You are attempting to access a block that is not stored"
  a.blocks[p]
end

