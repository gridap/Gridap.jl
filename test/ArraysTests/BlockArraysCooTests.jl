module BlockArraysCooTests

using Test
using Gridap.Arrays
using BlockArrays

blocks = [ [1 2; 3 4], [5 6; 7 8; 9 10] ]
blockids = [(1,1),(2,1)]
ax = (blockedrange([2,3]), blockedrange([2,4]))

a = BlockArrayCoo(blocks,blockids,ax)
@test a[Block(1),Block(1)] === blocks[1]
@test a[Block(1,1)] === blocks[1]
@test a[Block(1,2)] === a.zero_blocks[1]
@test a[BlockIndex((1,1),(2,1))] === blocks[1][2,1]
@test a[BlockIndex(1,2),BlockIndex(1,1)] === blocks[1][2,1]
@test a[2,1] === blocks[1][2,1]
@test a[3,2] === blocks[2][1,2]
@test axes(a) === ax
@test size(a) == (5,6)

@test is_zero_block(a,2,2) == true
@test is_zero_block(a,2,1) == false
@test is_nonzero_block(a,2,1) == true
@test is_zero_block(a,Block(1,2)) == true
@test is_zero_block(a,Block(1),Block(2)) == true

for (i,b) in enumerateblocks(a)
  @test a[i] === b
end

b21 = zeros(3,2)
getblock!(b21,a,Block(2,1))
@test b21 == a[Block(2,1)]

b12 = ones(2,4)
getblock!(b12,a,Block(1,2))
@test b12 == a[Block(1,2)]

blocks = [ [1,2,3] ]
blockids = [(2,)]
ax = (blockedrange([2,3,4]),)
a = BlockArrayCoo(blocks,blockids,ax)

@test a[Block(1)] === a.zero_blocks[1]
@test a[Block(2)] === blocks[1]
@test a[Block(3)] === a.zero_blocks[2]
@test a[BlockIndex(2,3)] === blocks[1][3]

using LinearAlgebra

const BlockMatrixCoo = BlockArrayCoo{T,2} where T
const BlockVectorCoo = BlockArrayCoo{T,1} where T

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

@static if VERSION >= v"1.3"

  function LinearAlgebra.mul!(c::BlockVectorCoo,a::BlockMatrixCoo,b::BlockVectorCoo,α,β)
    for I in 1:blocksize(a,1)
      cI = c[Block(I)]
      for i in eachindex(cI)
        cI[i] = β*cI[i]
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

  function LinearAlgebra.mul!(c::BlockMatrixCoo,a::BlockMatrixCoo,b::BlockMatrixCoo,α,β)
    for I in 1:blocksize(a,1)
      for J in 1:blocksize(b,2)
        cIJ = c[Block(I,J)]
        for ij in eachindex(cIJ)
          cIJ[ij] = β*cIJ[ij]
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

else

  function LinearAlgebra.mul!(c::BlockVectorCoo,a::BlockMatrixCoo,b::BlockVectorCoo)
    fill!(c,zero(eltype(c)))
    for I in 1:blocksize(a,1)
      for J in 1:blocksize(a,2)
        if is_nonzero_block(a,Block(I,J)) && is_nonzero_block(b,Block(J))
          aIJ = a[Block(I,J)]
          bJ = b[Block(J)]
          cI = c[Block(I)]
          for j in 1:size(aIJ,2)
            for i in 1:size(aIJ,1)
              cI[i] += aIJ[i,j]*bJ[j]
            end
          end
        end
      end
    end
    c
  end

  function LinearAlgebra.mul!(c::BlockMatrixCoo,a::BlockMatrixCoo,b::BlockMatrixCoo)
    fill!(c,zero(eltype(c)))
    for I in 1:blocksize(a,1)
      for J in 1:blocksize(b,2)
        for K in 1:blocksize(a,2)
          if is_nonzero_block(a,Block(I,K)) && is_nonzero_block(b,Block(K,J))
            aIK = a[Block(I,K)]
            bKJ = b[Block(K,J)]
            cIJ = c[Block(I,J)]
            for i in 1:size(aIK,1)
              for j in 1:size(bKJ,2)
                for k in 1:size(aIK,2)
                  cIJ[i,j] += aIK[i,k]*bKJ[k,j]
                end
              end
            end
          end
        end
      end
    end
    c
  end

end

function Base.fill!(a::BlockArrayCoo,v)
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

blocks = [ [1 2; 3 4], [5 6 7 8; 9 10 11 12; 13 14 15 16], [1 2 3 4; 5 6 7 8], [1 2 3; 4 5 6; 7 8 9] ]
blockids = [(1,1),(2,2),(1,2),(3,3)]
ax = (blockedrange([2,3,3]), blockedrange([2,4,3]))
a = BlockArrayCoo(blocks,blockids,ax)

blocks = [ 10*[1,2], 20*[1,2,3] ]
blockids = [(1,),(3,)]
ax = (blockedrange([2,4,3]),)
b = BlockArrayCoo(blocks,blockids,ax)

c = a*b
@test axes(c,1) === axes(a,1)
@test blocksize(c) == (3,)
@test Array(a)*Array(b) == c

mul!(c,a,b)
@test axes(c,1) === axes(a,1)
@test blocksize(c) == (3,)
@test Array(a)*Array(b) == c

b = Transpose(a)
c = a*b
@test axes(c,1) === axes(a,1)
@test axes(c,2) === axes(a,1)
@test blocksize(c) == (3,3)
@test Array(a)*Array(b) == c

mul!(c,a,b)
@test axes(c,1) === axes(a,1)
@test axes(c,2) === axes(a,1)
@test blocksize(c) == (3,3)
@test Array(a)*Array(b) == c

end # module
