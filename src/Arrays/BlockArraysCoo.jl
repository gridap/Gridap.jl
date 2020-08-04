
struct BlockArrayCoo{T,N,A,X} <: AbstractBlockArray{T,N}
  blocks::Vector{A}
  blockids::Vector{NTuple{N,Int}}
  axes::X
  ptrs::Array{Int,N}
  zero_blocks::Vector{A}
  function BlockArrayCoo(
    blocks::Vector{A},
    blockids::Vector{NTuple{N,Int}},
    axes::NTuple{N},
    ptrs::Array{Int,N}=_compute_ptrs(blockids,axes),
    zero_blocks::Vector{A}=_compute_zero_blocks(A,ptrs,axes)) where {T,N,A<:AbstractArray{T,N}}

    X = typeof(axes)
    new{T,N,A,X}(blocks,blockids,axes,ptrs,zero_blocks)
  end
end

const BlockMatrixCoo = BlockArrayCoo{T,2} where T
const BlockVectorCoo = BlockArrayCoo{T,1} where T

function _compute_ptrs(blockids,axes)
  s = map(i->first(blocksize(i)),axes)
  ptrs = zeros(Int,s)
  for (i,c) in enumerate(blockids)
    ptrs[c...] = i
  end
  m = 1
  for (k,p) in enumerate(ptrs)
    if p ==0
      ptrs[k] = -m
      m += 1
    end
  end
  ptrs
end

function _compute_zero_blocks(::Type{A},ptrs,axes) where A
  cis = CartesianIndices(ptrs)
  zero_blocks = A[]
  for ci in cis
    p = ptrs[ci]
    if p<0
      i = Tuple(ci)
      s = map((k,l)->length(k[Block(l)]),axes,i)
      block = A(undef,s)
      fill!(block,zero(eltype(A)))
      push!(zero_blocks,block)
    end
  end
  zero_blocks
end

function BlockArrays.getblock(a::BlockArrayCoo,i::Integer...)
  p = a.ptrs[i...]
  if p>0
    a.blocks[p]
  else
    a.zero_blocks[-p]
  end
end

@inline is_nonzero_block(a,b...) = ! is_zero_block(a,b...)

function is_zero_block(a,b::Block)
  i = convert(Tuple,b)
  is_zero_block(a,i...)
end

function is_zero_block(a,b::Block...)
  i = map(Int,b)
  is_zero_block(a,i...)
end

function is_zero_block(a,b::CartesianIndex)
  i = Tuple(b)
  is_zero_block(a,i...)
end

function is_zero_block(a::BlockArrayCoo,i::Integer...)
  p = a.ptrs[i...]
  @assert p != 0
  p < 0
end

function BlockArrays.getblock!(block,a::BlockArrayCoo,i::Integer...)
  p = a.ptrs[i...]
  if p>0
    block .= a.blocks[p]
  else
    block .= a.zero_blocks[-p]
  end
end

function BlockArrays.eachblock(a::BlockArrayCoo)
  cis = CartesianIndices(blocksize(a))
  blocks = map(ci->Block(Tuple(ci)),cis)
  ( a[block] for block in blocks  )
end

function enumerateblocks(a)
  cis = CartesianIndices(blocksize(a))
  blocks = map(ci->Block(Tuple(ci)),cis)
  zip(blocks,eachblock(a))
end

Base.size(a::BlockArrayCoo) = map(length,Base.axes(a))

Base.axes(a::BlockArrayCoo) = a.axes

function Base.getindex(a::BlockArrayCoo,i::Integer...)
  s = map(findblockindex,a.axes,i)
  ai = a[s...]
  ai
end

function Base.:*(a::BlockArrayCoo,b::BlockArrayCoo)
  c = _mul_block_result(a,b)
  mul!(c,a,b)
end

function _mul_block_result(a::BlockMatrixCoo,b::BlockVectorCoo)
  @assert blocksize(a,2) == blocksize(b,1)
  T = promote_type(eltype(a),eltype(b))
  blocks = Vector{T}[]
  blockids = Tuple{Int}[]
  for i in 1:blocksize(a,1)
    for j in 1:blocksize(a,2)
      if is_nonzero_block(a,Block(i,j)) && is_nonzero_block(b,Block(j))
        block = zeros(T,size(a[Block(i,j)],1))
        push!(blocks,block)
        push!(blockids,(i,))
        break
      end
    end
  end
  axs = (axes(a)[1],)
  BlockArrayCoo(blocks,blockids,axs)
end

function _mul_block_result(a::BlockMatrixCoo,b::BlockMatrixCoo)
  @assert blocksize(a,2) == blocksize(b,1)
  T = promote_type(eltype(a),eltype(b))
  blocks = Matrix{T}[]
  blockids = Tuple{Int,Int}[]
  for i in 1:blocksize(a,1)
    for j in 1:blocksize(b,2)
      for k in 1:blocksize(a,2)
        if is_nonzero_block(a,Block(i,k)) && is_nonzero_block(b,Block(k,j))
          block = zeros(T,size(a[Block(i,k)],1),size(b[Block(k,j)],2))
          push!(blocks,block)
          push!(blockids,(i,j))
          break
        end
      end
    end
  end
  axs = (axes(a)[1],axes(b)[2])
  BlockArrayCoo(blocks,blockids,axs)
end

#@static if VERSION >= v"1.3"

function LinearAlgebra.mul!(c::BlockVectorCoo,a::BlockMatrixCoo,b::BlockVectorCoo)
  mul!(c,a,b,1,0)
end

function LinearAlgebra.mul!(c::BlockMatrixCoo,a::BlockMatrixCoo,b::BlockMatrixCoo)
  mul!(c,a,b,1,0)
end

  function LinearAlgebra.mul!(c::BlockVectorCoo,a::BlockMatrixCoo,b::BlockVectorCoo,α::Number,β::Number)
    for I in 1:blocksize(a,1)
      cI = c[Block(I)]
      if is_nonzero_block(c,Block(I))
        for i in eachindex(cI)
          cI[i] = β*cI[i]
        end
      end
      for J in 1:blocksize(a,2)
        if is_nonzero_block(a,Block(I,J)) && is_nonzero_block(b,Block(J))
          aIJ = a[Block(I,J)]
          bJ = b[Block(J)]
          mul!(cI,aIJ,bJ,α,1)
        end
      end
    end
    c
  end

  function LinearAlgebra.mul!(c::BlockMatrixCoo,a::BlockMatrixCoo,b::BlockMatrixCoo,α::Number,β::Number)
    for I in 1:blocksize(a,1)
      for J in 1:blocksize(b,2)
        cIJ = c[Block(I,J)]
        if is_nonzero_block(c,Block(I,J))
          for ij in eachindex(cIJ)
            cIJ[ij] = β*cIJ[ij]
          end
        end
        for K in 1:blocksize(a,2)
          if is_nonzero_block(a,Block(I,K)) && is_nonzero_block(b,Block(K,J))
            aIK = a[Block(I,K)]
            bKJ = b[Block(K,J)]
            mul!(cIJ,aIK,bKJ,α,1)
          end
        end
      end
    end
    c
  end

#else
#
#  function LinearAlgebra.mul!(c::BlockVectorCoo,a::BlockMatrixCoo,b::BlockVectorCoo)
#    fill!(c,zero(eltype(c)))
#    for I in 1:blocksize(a,1)
#      for J in 1:blocksize(a,2)
#        if is_nonzero_block(a,Block(I,J)) && is_nonzero_block(b,Block(J))
#          aIJ = a[Block(I,J)]
#          bJ = b[Block(J)]
#          cI = c[Block(I)]
#          for j in 1:size(aIJ,2)
#            for i in 1:size(aIJ,1)
#              cI[i] += aIJ[i,j]*bJ[j]
#            end
#          end
#        end
#      end
#    end
#    c
#  end
#
#  function LinearAlgebra.mul!(c::BlockMatrixCoo,a::BlockMatrixCoo,b::BlockMatrixCoo)
#    fill!(c,zero(eltype(c)))
#    for I in 1:blocksize(a,1)
#      for J in 1:blocksize(b,2)
#        for K in 1:blocksize(a,2)
#          if is_nonzero_block(a,Block(I,K)) && is_nonzero_block(b,Block(K,J))
#            aIK = a[Block(I,K)]
#            bKJ = b[Block(K,J)]
#            cIJ = c[Block(I,J)]
#            for i in 1:size(aIK,1)
#              for j in 1:size(bKJ,2)
#                for k in 1:size(aIK,2)
#                  cIJ[i,j] += aIK[i,k]*bKJ[k,j]
#                end
#              end
#            end
#          end
#        end
#      end
#    end
#    c
#  end
#
#end

function Base.fill!(a::BlockArrayCoo,v)
  @notimplementedif v != zero(v)
  for b in a.blocks
    fill!(b,v)
  end
  a
end

function LinearAlgebra.Transpose(a::BlockMatrixCoo)
  blocks = [ Transpose(block) for block in a.blocks ]
  zero_blocks = [ Transpose(block) for block in a.zero_blocks ]
  blockids = [ (j,i) for (i,j) in a.blockids ]
  ax,ay = axes(a)
  axs = (ay,ax)
  ptrs = collect(Transpose(a.ptrs))
  BlockArrayCoo(blocks,blockids,axs,ptrs,zero_blocks)
end

function Base.copy(a::BlockArrayCoo)
  blocks = copy(a.blocks)
  blockids = copy(a.blockids)
  axes = copy.(a.axes)
  ptrs = copy(a.ptrs)
  zero_blocks = copy(a.zero_blocks)
  BlockArrayCoo(blocks,blockids,axes,ptrs,zero_blocks)
end

function Base.copyto!(a::BlockArrayCoo,b::BlockArrayCoo)
  for p in 1:length(a.blocks)
    copyto!(a.blocks[p],b.blocks[p])
  end
 a
end
