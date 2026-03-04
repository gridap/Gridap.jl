
# Make Base.:* behave like a Map

return_value(::typeof(*),a::Number,b::AbstractArray{<:Number}) = a*b

function return_cache(::typeof(*),a::Number,b::AbstractArray{<:Number})
  T = typeof(a*testitem(b))
  N = ndims(b)
  CachedArray(T,N)
end

function evaluate!(cache,::typeof(*),a::Number,b::AbstractArray{<:Number})
  setsize!(cache,size(b))
  c = cache.array
  mul!(c,b,a)
  return c
end

return_value(::typeof(*),a::AbstractArray{<:Number},b::Number) = a*b
return_cache(::typeof(*),a::AbstractArray{<:Number},b::Number) = return_cache(*,b,a)
evaluate!(cache,::typeof(*),a::AbstractArray{<:Number},b::Number) = evaluate!(cache,*,b,a)

function return_cache(::typeof(*),a::AbstractArray{<:Number},b::AbstractArray{<:Number})
  T = typeof(testitem(a)*testitem(b))
  N = ifelse(isa(b,AbstractVector),1,2)
  CachedArray(T,N)
end

function evaluate!(cache,::typeof(*),a::AbstractArray{<:Number},b::AbstractArray{<:Number})
  setsize_op!(*,cache,a,b)
  c = cache.array
  mul!(c,a,b)
  return c
end

function return_value(::typeof(*),a::ArrayBlock{A,2},b::ArrayBlock{B,1}) where {A,B}
  ri = return_value(*,testvalue(A),testvalue(B))
  array = Vector{typeof(ri)}(undef,size(a.array,1))
  touched = fill(false,size(a.array,1))
  ArrayBlock(array,touched)
end

function return_cache(::typeof(*),a::ArrayBlock,b::ArrayBlock)
  c1 = CachedArray(a*b)
  c2 = return_cache(unwrap_cached_array,c1)
  return c1, c2
end

function evaluate!(cache,::typeof(*),a::ArrayBlock,b::ArrayBlock)
  c1, c2 = cache
  setsize_op!(*,c1,a,b)
  c = evaluate!(c2,unwrap_cached_array,c1)
  mul!(c,a,b)
  c
end

# MulAddMap: Cached version of `mul!(d,a,b,α,β)`

struct MulAddMap{T} <: Map
  α::T
  β::T
end

function return_cache(k::MulAddMap,a,b,c)
  d = a*b+c
  CachedArray(d)
end

function evaluate!(cache,k::MulAddMap,a,b,c)
  setsize!(cache,size(c))
  d = cache.array
  copyto!(d,c)
  iszero(k.α) && isone(k.β) && return d
  mul!(d,a,b,k.α,k.β)
  d
end

function return_value(k::MulAddMap,a::ArrayBlock,b::ArrayBlock,c::ArrayBlock)
  x = return_value(*,a,b)
  return_value(+,x,c)
end

function return_cache(k::MulAddMap,a::ArrayBlock,b::ArrayBlock,c::ArrayBlock)
  c1 = CachedArray(a*b+c)
  c2 = return_cache(unwrap_cached_array,c1)
  (c1,c2)
end

function evaluate!(cache,k::MulAddMap,a::ArrayBlock,b::ArrayBlock,c::ArrayBlock)
  c1,c2 = cache
  setsize_op!(copy,c1,c)
  setsize_op!(*,c1,a,b)
  d = evaluate!(c2,unwrap_cached_array,c1)
  copyto!(d,c)
  iszero(k.α) && isone(k.β) && return d
  mul!(d,a,b,k.α,k.β)
  d
end

# Assembly Maps: AddEntriesMap and TouchEntriesMap

struct AddEntriesMap{F} <: Map
  combine::F
end

function evaluate!(cache,k::AddEntriesMap,A,v,i,j)
  add_entries!(k.combine,A,v,i,j)
end

function evaluate!(cache,k::AddEntriesMap,A,v,i)
  add_entries!(k.combine,A,v,i)
end

struct TouchEntriesMap <: Map end

function evaluate!(cache,k::TouchEntriesMap,A,v,i,j)
  add_entries!(+,A,nothing,i,j)
end

function evaluate!(cache,k::TouchEntriesMap,A,v,i)
  add_entries!(+,A,nothing,i)
end

for T in (:AddEntriesMap,:TouchEntriesMap)
  @eval begin

    function return_cache(k::$T,A,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
      qs = findall(v.touched)
      i, j = Tuple(first(qs))
      cij = return_cache(k,A,v.array[i,j],I.array[i],J.array[j])
      ni,nj = size(v.touched)
      cache = Matrix{typeof(cij)}(undef,ni,nj)
      for j in 1:nj
        for i in 1:ni
          if v.touched[i,j]
            cache[i,j] = return_cache(k,A,v.array[i,j],I.array[i],J.array[j])
          end
        end
      end
      cache
    end

    function evaluate!(cache, k::$T,A,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
      ni,nj = size(v.touched)
      for j in 1:nj
        for i in 1:ni
          if v.touched[i,j]
            evaluate!(cache[i,j],k,A,v.array[i,j],I.array[i],J.array[j])
          end
        end
      end
    end

    function return_cache(k::$T,A,v::VectorBlock,I::VectorBlock)
      qs = findall(v.touched)
      i = first(qs)
      ci = return_cache(k,A,v.array[i],I.array[i])
      ni = length(v.touched)
      cache = Vector{typeof(ci)}(undef,ni)
      for i in 1:ni
        if v.touched[i]
          cache[i] = return_cache(k,A,v.array[i],I.array[i])
        end
      end
      cache
    end

    function evaluate!(cache, k::$T,A,v::VectorBlock,I::VectorBlock)
      ni = length(v.touched)
      for i in 1:ni
        if v.touched[i]
          evaluate!(cache[i],k,A,v.array[i],I.array[i])
        end
      end
    end
  end

  for MT in (:MatrixBlock,:MatrixBlockView)
    Aij = (MT == :MatrixBlock) ? :(A.array[i,j]) : :(A[i,j])
    @eval begin

      function return_cache(k::$T,A::$MT,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
        qs = findall(v.touched)
        i, j = Tuple(first(qs))
        cij = return_cache(k,$Aij,v.array[i,j],I.array[i],J.array[j])
        ni,nj = size(v.touched)
        cache = Matrix{typeof(cij)}(undef,ni,nj)
        for j in 1:nj
          for i in 1:ni
            if v.touched[i,j]
              cache[i,j] = return_cache(k,$Aij,v.array[i,j],I.array[i],J.array[j])
            end
          end
        end
        cache
      end

      function evaluate!(cache,k::$T,A::$MT,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
        ni,nj = size(v.touched)
        for j in 1:nj
          for i in 1:ni
            if v.touched[i,j]
              evaluate!(cache[i,j],k,$Aij,v.array[i,j],I.array[i],J.array[j])
            end
          end
        end
      end

    end # @eval
  end # for MT

  for VT in (:VectorBlock,:VectorBlockView)
    Ai = (VT == :VectorBlock) ? :(A.array[i]) : :(A[i])
    @eval begin

      function return_cache(k::$T,A::$VT,v::VectorBlock,I::VectorBlock)
        qs = findall(v.touched)
        i = first(qs)
        ci = return_cache(k,$Ai,v.array[i],I.array[i])
        ni = length(v.touched)
        cache = Vector{typeof(ci)}(undef,ni)
        for i in 1:ni
          if v.touched[i]
            cache[i] = return_cache(k,$Ai,v.array[i],I.array[i])
          end
        end
        cache
      end

      function evaluate!(cache, k::$T,A::$VT,v::VectorBlock,I::VectorBlock)
        ni = length(v.touched)
        for i in 1:ni
          if v.touched[i]
            evaluate!(cache[i],k,$Ai,v.array[i],I.array[i])
          end
        end
      end

    end # @eval
  end # for VT
end # end for T

###########################################################################################
# Local solver maps

"""
    LocalSolveMap(; factorize! = lu!, pivot = RowMaximum())

A map for solving local linear systems, relying on a factorization method.

Given a left-hand-side matrix `mat` and a set of N right-hand-side arrays `lhs`, 
returns an N-Tuple of arrays containing the solutions to the linear systems defined by 

Each system is given by `A*x_i = b_i`, and the solution is computed as 
`x_i = ldiv!(factorize!(A,pivot),b_i)`

"""
struct LocalSolveMap{A,B} <: Map
  factorize! :: A
  pivot :: B
end

LocalSolveMap(; factorize! = lu!, pivot = RowMaximum()) = LocalSolveMap(factorize!, pivot)

Arrays.return_value(k::LocalSolveMap, matvec::Tuple) = return_value(k,matvec[1],matvec[2])
Arrays.return_cache(k::LocalSolveMap, matvec::Tuple) = return_cache(k,matvec[1],matvec[2])
Arrays.evaluate!(cache, k::LocalSolveMap, matvec::Tuple) = evaluate!(cache,k,matvec[1],matvec[2])

Arrays.return_value(k::LocalSolveMap, mat, vec) = similar(vec)

function Arrays.return_cache(k::LocalSolveMap, mat, vec)
  CachedArray(mat), CachedArray(vec)
end

function Arrays.evaluate!(cache,k::LocalSolveMap, mat, vec)
  cmat, cvec = cache

  setsize!(cmat, size(mat))
  copy!(cmat.array, mat)
  f = k.factorize!(cmat.array,k.pivot;check=false)
  @check issuccess(f) "Factorization failed"

  setsize!(cvec, size(vec))
  ldiv!(cvec.array,f,vec)

  return cvec.array
end

function Arrays.return_value(k::LocalSolveMap, mat::MatrixBlock, vec::MatrixBlock)
  ntuple(i -> similar(get_array(vec)[i]), size(vec,2))
end

function Arrays.return_cache(k::LocalSolveMap, mat::MatrixBlock, vec::MatrixBlock)
  @check size(mat,1) == size(mat,2) == size(vec,1) == 1
  CachedArray(get_array(mat)[1,1]), ntuple(i -> CachedArray(get_array(vec)[i]), size(vec,2))
end

function Arrays.evaluate!(cache,k::LocalSolveMap, mat::MatrixBlock, vec::MatrixBlock)
  @check size(mat,1) == size(mat,2) == size(vec,1) == 1
  cmat, cvec = cache

  m = get_array(mat)[1,1]
  setsize!(cmat, size(m))
  copy!(cmat.array, m)
  f = k.factorize!(cmat.array,k.pivot;check=false)
  @check issuccess(f) "Factorization failed"

  v = get_array(vec)
  for i in eachindex(v)
    setsize!(cvec[i], size(v[i]))
    ldiv!(cvec[i].array,f,v[i])
  end

  return ntuple(i -> cvec[i].array, size(vec,2))
end

###########################################################################################

"""
    LocalPenaltySolveMap(; factorize! = lu!, pivot = RowMaximum())

A map for solving local constrained linear systems, relying on a factorization method.

Given a left-hand-side 2x2 block matrix matrix`mat` and a set of 2xN right-hand-side arrays `lhs`, 
returns an N-Tuple of arrays containing the solutions to the linear systems. 

Each system is given by `A*[x_i; λ_i] = b_i`, where `A = [App, Aλp; Apλ, 0]` is the
augmented matrix, and `b_i = [Bp; Bλ]` is the right-hand side vector. The solution is 
computed using a penalty method, as `x_i = ldiv!(factorize!(C,pivot),d_i)` with 
`C = App + μT * Apλ * Aλp` and `d_i = Bp + μT * Apλ * Bλ`, where `μT` is a penalty parameter.
The penalty parameter μT is heuristically chosen as `μT = norm(App)/norm(Apλ*Aλp)`.

"""
struct LocalPenaltySolveMap{A,B} <: Map
  factorize! :: A
  pivot :: B
end

LocalPenaltySolveMap(; factorize! = lu!, pivot = RowMaximum()) = LocalPenaltySolveMap(factorize!, pivot)

function Arrays.evaluate!(cache::Nothing, k::LocalPenaltySolveMap, mats::Tuple)
  lhs, rhs = mats
  evaluate!(cache, k, lhs, rhs)
end

function Arrays.evaluate!(cache::Nothing, k::LocalPenaltySolveMap, lhs, rhs::VectorBlock)
  @check size(lhs) == (2,2)
  @check size(rhs) == (2,)

  App, Aλp, Apλ, _ = get_array(lhs)
  Bp, Bλ = get_array(rhs)

  # μT = norm(App)/norm(Apλ*Aλp) is a heuristic choice for the penalty parameter
  if isone(size(Apλ,2))
    μT = tr(App)/norm(Apλ)^2 # Single constraint
  else
    μT = tr(App)/norm(Apλ*Aλp) # Multiple constraints
  end
  
  # App = App + μT * Apλ * Aλp
  mul!(App, Apλ, Aλp, μT, 1)
  
  # Bp = Bp + μT * Apλ * Bλ
  mul!(Bp, Apλ, Bλ, μT, 1)

  Ainv = k.factorize!(App,k.pivot;check=false)
  @check issuccess(Ainv) "Factorization failed"

  Ru = ldiv!(Ainv,Bp)
  return (Ru,)
end

function Arrays.evaluate!(cache::Nothing, k::LocalPenaltySolveMap, lhs, rhs)
  @check size(lhs) == (2,2)
  @check size(rhs) == (2,2)

  App, Aλp, Apλ, _ = get_array(lhs)
  BpΩ, BλΩ, BpΓ, _ = get_array(rhs)

  # μT = norm(App)/norm(Apλ*Aλp) is a heuristic choice for the penalty parameter
  if isone(size(Apλ,1))
    μT = tr(App)/norm(Apλ)^2 # Single constraint
  else
    μT = tr(App)/norm(Apλ*Aλp) # Multiple constraints
  end
  
  # App = App + μT * Apλ * Aλp
  mul!(App, Apλ, Aλp, μT, 1)
  
  # BpΩ = BpΩ + μT * Apλ * BλΩ
  mul!(BpΩ, Apλ, BλΩ, μT, 1)

  Ainv = k.factorize!(App,k.pivot;check=false)
  @check issuccess(Ainv) "Factorization failed"

  RuΩ = ldiv!(Ainv,BpΩ)
  RuΓ = ldiv!(Ainv,BpΓ)

  return RuΩ, RuΓ
end
