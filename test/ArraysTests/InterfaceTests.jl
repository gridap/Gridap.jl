module InterfaceTests

using Test
using Gridap.Arrays

import Gridap.Arrays: getindex!

a = rand(20,12)

test_array(a,a)
test_array(a,a,≈)


# LinearIndex index style
#1D array
struct IndexlinVector{T} <: AbstractVector{T}
  v::Vector{T}
end
Base.size(lv::IndexlinVector) = size(lv.v)
Base.IndexStyle(::Type{<:IndexlinVector}) = IndexLinear()
Base.getindex(lv::IndexlinVector, i::Integer) = lv.v[i]

v = [1, 2, 3, 4]
lv = IndexlinVector(v)
cache = array_cache(lv)

# default indexing without cache
@test Arrays.uses_hash(typeof(lv)) == Val(false)
@test isnothing(cache)
@test_nowarn getindex!(nothing, lv, CartesianIndex(3))
@test_nowarn getindex!(nothing, lv, 3)

# indexing using cache
Arrays.getindex!(cache, lv::IndexlinVector, i::Integer) = begin @warn("passed"); lv.v[i] end
@test_warn "passed" Arrays.getindex!(nothing, lv, CartesianIndex(3))
@test_warn "passed" Arrays.getindex!(nothing, lv, 3)
@test_nowarn Arrays.getindex!(nothing, lv, 1:2) # cache skipped
test_array(v, lv)


#2D array
struct IndexlinMatrix{T} <: AbstractMatrix{T}
  m::Matrix{T}
end
Base.size(lm::IndexlinMatrix) = size(lm.m)
Base.IndexStyle(::Type{<:IndexlinMatrix}) = IndexLinear()
Base.getindex(lm::IndexlinMatrix, i::Integer) = lm.m[i]

m = [1 2; 3 4]
lm = IndexlinMatrix(m)
cache = array_cache(lm)

# default indexing without cache
@test Arrays.uses_hash(typeof(lm)) == Val(false)
@test isnothing(cache)
@test_nowarn getindex!(nothing, lm, 1, 2)
@test_nowarn getindex!(nothing, lm, CartesianIndex(1, 2))
@test_nowarn getindex!(nothing, lm, 3)

# indexing using cache
Arrays.getindex!(cache, lm::IndexlinMatrix, i::Integer) = begin @warn("passed"); lm.m[i] end
@test_warn "passed" Arrays.getindex!(nothing, lm, 3)
@test_warn "passed" Arrays.getindex!(nothing, lm, 1, 2)
@test_warn "passed" Arrays.getindex!(nothing, lm, CartesianIndex(1, 2))
@test_warn "passed" Arrays.getindex!(nothing, lm, CartesianIndex(1), CartesianIndex(2))
@test_nowarn Arrays.getindex!(nothing, lm, 1:2,2) # cache skipped
test_array(m, lm)


# CartesianIndex index style
#1D array
struct IndexcartVector{T} <: AbstractVector{T}
  v::Vector{T}
end
Base.size(cv::IndexcartVector) = size(cv.v)
Base.IndexStyle(::Type{<:IndexcartVector}) = IndexCartesian()
Base.getindex(cv::IndexcartVector, i::Integer) = cv.v[i]

v = [1, 2, 3, 4]
cv = IndexcartVector(v)
cache = array_cache(cv)

# default indexing without cache
@test Arrays.uses_hash(typeof(cv)) == Val(false)
@test isnothing(cache)
@test_nowarn getindex!(nothing, cv, 3)
@test_nowarn getindex!(nothing, cv, CartesianIndex(3))

# indexing using cache
Arrays.getindex!(cache, cv::IndexcartVector, i::Integer) = begin @warn("passed"); cv.v[i] end
@test_warn "passed" Arrays.getindex!(nothing, cv, 3)
@test_warn "passed" Arrays.getindex!(nothing, cv, CartesianIndex(3))
@test_nowarn Arrays.getindex!(nothing, cv, 1:2) # cache skipped
test_array(v, cv)


#2D array
struct IndexcartMatrix{T} <: AbstractMatrix{T}
  m::Matrix{T}
end
Base.size(cm::IndexcartMatrix) = size(cm.m)
Base.IndexStyle(::Type{<:IndexcartMatrix}) = IndexCartesian()
Base.getindex(cm::IndexcartMatrix, i::Integer, j::Integer) = cm.m[i,j]

m = [1 2; 3 4]
cm = IndexcartMatrix(m)
cache = array_cache(cm)

# default indexing without cache
@test Arrays.uses_hash(typeof(cm)) == Val(false)
@test isnothing(cache)
@test_nowarn getindex!(nothing, cm, 1, 2)
@test_nowarn getindex!(nothing, cm, CartesianIndex(1, 2))
@test_nowarn getindex!(nothing, cm, 3)

# indexing using cache
Arrays.getindex!(cache, cm::IndexcartMatrix, i::Integer, j::Integer) = begin @warn("passed"); cm.m[i,j] end
@test_warn "passed" Arrays.getindex!(nothing, cm, 1, 2)
@test_warn "passed" Arrays.getindex!(nothing, cm, CartesianIndex(1, 2))
@test_warn "passed" Arrays.getindex!(nothing, cm, CartesianIndex(1), CartesianIndex(2))
@test_warn "passed" Arrays.getindex!(nothing, cm, 3)
@test_nowarn Arrays.getindex!(nothing, cm, 1:2,2) # cache skipped
test_array(m, cm)


t = (1,2.0,1im,[3,4])
@test testvalue(t) == (0,0.0,0im,Int[])

for n in 0:10
  t = ntuple(i->rand(),n)
  @test testvalue(t) == ntuple(i->0.0,n)
end

end # module
