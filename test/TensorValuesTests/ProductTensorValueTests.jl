module ProductTensorValueTests

using Gridap.TensorValues
using Gridap.TensorValues: _minimal_tensortype_with_shape, _factor_eltype, stored_comp_indices
using Gridap.TensorValues: _prod_counts, _prod_eltype, _rm_size_variance, _Size_to_size
using Gridap.TensorValues: get_indep_components, num_stored_components
using Gridap.TensorValues: _indep_comp_cartesian_indices, is_full_rank
using Gridap.TensorValues: _stored_comp_cartesian_indices
using Gridap.TensorValues: stored_comp_indices, _stored_comp_weights, _dense_flat_inner
using Gridap.TensorValues: full_rank, is_opposit_variance
using Gridap.Arrays: get_array
using Test
using StaticArrays

# Factor types, given as the eltype-free tags the product stores. Kept small:
# the conversions to `SArray` unroll `prod(S)` components, so a pair of
# high-order factors costs a lot to compile.
tags = [
  VectorValue{2}, VectorValue{3},
  CoVectorValue{cov2}, ContraVectorValue{con3},
  TensorValue{2,3},
  SymTensorValue{2}, SymTensorValue{3},
  SymTracelessTensorValue{2}, SymTracelessTensorValue{3},
  SkewSymTensorValue{3},
  ExteriorFormValue{2,cov3},
  NTensorValue{Tuple{2,2,2}},
  SymNTensorValue{Tuple{cov2,cov2,cov2}}
 ]

# Canonicalisation of the factor tags

@test _minimal_tensortype_with_shape(VectorValue{3,Float64}) === VectorValue{3}
@test _minimal_tensortype_with_shape(SymTensorValue{3,Float64,6}) === SymTensorValue{3}
@test _minimal_tensortype_with_shape(SymTensorValue{3}) === SymTensorValue{3}
@test _minimal_tensortype_with_shape(ThirdOrderTensorValue{2,3,4,Int,24}) === ThirdOrderTensorValue{2,3,4}
@test _minimal_tensortype_with_shape(SymNTensorValue{Tuple{cov2,cov2,cov2},Int,3,4}) === SymNTensorValue{Tuple{cov2,cov2,cov2}}
@test _minimal_tensortype_with_shape(CoVectorValue{cov3,Float64,3}) === ExteriorFormValue{1,cov3}

@test _factor_eltype(VectorValue{3}) === nothing
@test _factor_eltype(ThirdOrderTensorValue{2,3,4}) === nothing
@test _factor_eltype(CoVectorValue{cov3}) === nothing
@test _factor_eltype(VectorValue{3,Float64}) === Float64
@test _factor_eltype(SymTensorValue{3,Int}) === Int
@test _factor_eltype(CoVectorValue{cov3,Float64,3}) === Float64

# The counts and the component type are constant-folded from the factor types,
# which is what keeps `_get_data` free of the `divrem` and construction free of
# allocation; `Val` fails to infer as soon as one of them is only known at run
# time

_counts_val() = Val(_prod_counts(VectorValue{3}, SymTracelessTensorValue{3}))
_eltype_val() = Val(_prod_eltype(VectorValue{3,Float64}, VectorValue{2}))

@test @inferred(_counts_val()) === Val((3,6,3,5))
@test @inferred(_eltype_val()) === Val(Float64)

a, b = VectorValue(1.0,2.0), SymTensorValue(1,2,3)
V = ProductTensorValue{VectorValue{2},SymTensorValue{2}}
P = V{Float64,Tuple{2,2,2},3,6}
@test typeof(product_tensor(a,b)) === P

# A component type given in a factor is not kept in the factor tags, it gives
# the component type of the product instead

@test typeof(@inferred(ProductTensorValue{VectorValue{2,Float64},SymTensorValue{2}}(ntuple(Int,6)))) === P
@test typeof(ProductTensorValue{VectorValue{2},SymTensorValue{2,Float64,3}}(ntuple(Int,6))) === P
@test typeof(ProductTensorValue{VectorValue{2,Int},SymTensorValue{2,Float64,3}}(ntuple(Int,6))) === P
@test typeof(V(ntuple(Float64,6))) === P
@test typeof(V(ntuple(Int,6))) === change_eltype(P,Int)

# An explicitly given component type prevails over that of the factors, as long
# as the latter promotes to it

@test typeof(ProductTensorValue{VectorValue{2,Int},SymTensorValue{2,Int,3},Float64,Tuple{2,2,2},3,6}(ntuple(Int,6))) === P
@test_throws AssertionError ProductTensorValue{VectorValue{2,Float64},SymTensorValue{2},Int,Tuple{2,2,2},3,6}(ntuple(Int,6))

# `S`, `N` and `L` follow from the factor types, so the component type is the
# only other parameter to give

@test typeof(@inferred(V{Float64}(ntuple(Int,6)))) === P
@test typeof(V{Float64}(ntuple(Int,6)...)) === P
@test_throws AssertionError ProductTensorValue{VectorValue{2,Float64},SymTensorValue{2},Int}(ntuple(Int,6))

# `S`, `N` and `L` may still be spelled out, in part or in full, and are then
# checked against the ones the factors give

@test typeof(@inferred(V{Float64,Tuple{2,2,2}}(ntuple(Int,6)))) === P
@test typeof(V{Float64,Tuple{2,2,2}}(ntuple(Int,6)...)) === P
@test typeof(@inferred(V{Float64,Tuple{2,2,2},3}(ntuple(Int,6)))) === P
@test typeof(V{Float64,Tuple{2,2,2},3}(ntuple(Int,6)...)) === P

@test_throws AssertionError V{Float64,Tuple{2,2,3}}(ntuple(Int,6))
@test_throws AssertionError V{Float64,Tuple{2,2,2},2}(ntuple(Int,6))

# All the spellings build the same value

@test allequal((
  V(ntuple(Float64,6)),
  V{Float64}(ntuple(Int,6)),
  V{Float64}(ntuple(Int,6)...),
  V{Float64,Tuple{2,2,2}}(ntuple(Int,6)),
  V{Float64,Tuple{2,2,2}}(ntuple(Int,6)...),
  V{Float64,Tuple{2,2,2},3}(ntuple(Int,6)),
  V{Float64,Tuple{2,2,2},3}(ntuple(Int,6)...),
  V{Float64,Tuple{2,2,2},3,6}(ntuple(Int,6)),
  V{Float64,Tuple{2,2,2},3,6}(ntuple(Int,6)...)))

# Array constructors tests

full = [ Float64(i + 10j + 100k) for i in 1:2, j in 1:2, k in 1:2 ]
full[:,1,2] = full[:,2,1]

@test typeof(@inferred(V(full))) === P
@test get_array(V(full)) == full

@test allequal((
  V(full),
  V{Float64}(full),
  V{Float64,Tuple{2,2,2}}(full),
  V{Float64,Tuple{2,2,2},3}(full),
  V{Float64,Tuple{2,2,2},3,6}(full),
  V(SArray{Tuple{2,2,2}}(full)),
  V{Float64}(Int.(full))))

@test eltype(V(Int.(full))) === Int
@test eltype(ProductTensorValue{VectorValue{2,Float64},SymTensorValue{2}}(Int.(full))) === Float64
@test typeof(ProductTensorValue{VectorValue{0},SymTensorValue{2}}(zeros(0,2,2))) ===
      ProductTensorValue{VectorValue{0},SymTensorValue{2},Float64,Tuple{0,2,2},3,0}

@test_throws AssertionError V(zeros(2,2,3))
@test_throws AssertionError V(zeros(8))

# A full tensor of the right size but of offset indices would be read at the
# wrong components
@test_throws ArgumentError V(view(zeros(3,3,3), ntuple(_ -> Base.IdentityUnitRange(2:3), 3)...))

# The components that a factor makes dependent are never read, so a full tensor
# lacking the symmetry of its factor is read at its representatives only
@test Tuple(ProductTensorValue{SymTensorValue{2},VectorValue{1}}(reshape(Float64[1 2; 3 4],2,2,1))) ==
      Tuple(ProductTensorValue{SymTensorValue{2},VectorValue{1}}(1.,3.,4.))

# Storage layout: the independent components of the first factor vary fastest,
# so a product of two tensors without dependent components has the components
# of their `outer` product, in the same order.

p = product_tensor(VectorValue(1,2), VectorValue(10,20,30))
@test Tuple(p) == (10,20,20,40,30,60)
@test Tuple(p) == Tuple(VectorValue(1,2) ⊗ VectorValue(10,20,30))
@test size(p) == (2,3)
@test num_indep_components(p) == 6

p = product_tensor(VectorValue(1,2), SymTensorValue(10,20,30))
@test Tuple(p) == (10,20,20,40,30,60)
@test num_indep_components(p) == 6
@test num_components(p) == 8

# A factor that materialises its dependent components makes the product store
# them too, so that reading a component of the product stays a single read

q, v = SymTracelessTensorValue(1,2,3,4,5), VectorValue(1,10)
@test (num_indep_components(q), num_stored_components(typeof(q))) == (5,6)

# The stored count is `L`, and falls back to the independent count for every
# type that does not materialise dependent components

_stored_L(::Type{<:MultiValue{S,T,N,L}}) where {S,T,N,L} = L

@test num_stored_components(q) == num_stored_components(typeof(q))
@test num_stored_components(SymTracelessTensorValue{3}) == 6
@test num_stored_components(SymTensorValue{3,Float64}) == 6
@test num_stored_components(VectorValue{3}) == 3

# A traceless tensor in dimension zero stores nothing, its symmetric layout
# being empty too

@test num_stored_components(SymTracelessTensorValue{0}) == _stored_L(SymTracelessTensorValue{0,Int,0}) == 0
@test all(V -> num_stored_components(V) == _stored_L(make_concretetype(V{Int})), tags)
@test_throws "too abstract to count its components" num_stored_components(VectorValue)

pq = product_tensor(q,v)
@test num_stored_components(typeof(pq)) == _stored_L(typeof(pq)) == 12
@test num_indep_components(pq) == 10
@test length(Tuple(pq)) == 12
@test get_indep_components(pq) == ntuple(i -> indep_comp_getindex(pq,i), 10)
@test indep_comp_getindex(pq,7) == indep_comp_getindex(q,2)*indep_comp_getindex(v,2)
@test all(pq[i,j,k] == q[i,j]*v[k] for i in 1:3, j in 1:3, k in 1:2)
@test pq == outer(q,v)
@test typeof(@inferred(pq+pq)) === typeof(pq)
@test Tuple(pq+pq) == 2 .* Tuple(pq)
@test typeof(@inferred(2pq)) === typeof(pq)
@test iszero(pq - pq)

# The dependent component of a traceless tensor is stored last, in the layout of
# a symmetric tensor, which `stored_comp_indices` uses as its slot map

@test Tuple(q) == Tuple(SymTensorValue(1,2,3,4,5,-5))
@test stored_comp_indices(typeof(q)) === SymTensorValue{3,Int}(1,2,3,4,5,6)
@test stored_comp_indices(VectorValue{3,Int}) === VectorValue{3,Int}(1,2,3)

# A factor tag carrying only its size information is enough, so the slot map of
# a factor is read off the tag the product stores

@test stored_comp_indices(VectorValue{3}) === VectorValue{3,Int}(1,2,3)
@test stored_comp_indices(SymTracelessTensorValue{3}) === stored_comp_indices(typeof(q))
@test all(V -> stored_comp_indices(V) === stored_comp_indices(make_concretetype(V{Int})), tags)

# All the components of a 1×1 traceless tensor are dependent, so the product has
# no independent component to give its component type

p0 = product_tensor(zero(SymTracelessTensorValue{1,Float64}), v)
@test num_indep_components(p0) == 0
@test length(Tuple(p0)) == 2
@test eltype(p0) === Float64
@test iszero(p0)

# Type parameters and inference

a, b = VectorValue(1.0,2.0), SymTensorValue(1.0,2.0,3.0)
p = @inferred product_tensor(a,b)
@test typeof(p) === V{Float64,Tuple{2,2,2},3,6}
@test @inferred(zero(typeof(p))) == product_tensor(zero(a),b)
@test typeof(@inferred(p+p)) === typeof(p)
@test typeof(@inferred(2.0*p)) === typeof(p)
@test typeof(@inferred(2*p)) === typeof(p)
@test @inferred(p[1,2,2]) == a[1]*b[2,2]
@test typeof(@inferred(V{Float64,Tuple{2,2,2},3,6}(ntuple(Float64,6)))) === typeof(p)

@test change_eltype(p,Int) === V{Int,Tuple{2,2,2},3,6}
@test promote_type(typeof(p),change_eltype(p,Int)) === typeof(p)
@test convert(change_eltype(p,Float32),p) == p

# A product distributes the dimensions over its factors and changes the variance
# of each, keeping its factorisation
Pv = ProductTensorValue{CoVectorValue{cov2},SymTensorValue{con2}}
@test change_variance(Pv, Val((con2,cov2,cov2))) ===
      ProductTensorValue{ContraVectorValue{con2},SymTensorValue{cov2}}
@test change_variance(Pv{Float64}, Val((con2,cov2,cov2))) ===
      ProductTensorValue{ContraVectorValue{con2},SymTensorValue{cov2},Float64}
@test change_variance(Pv, inverse_variance(Tuple{cov2,con2,con2})) ===
      ProductTensorValue{ContraVectorValue{con2},SymTensorValue{cov2}}
# each factor applies its own contract to its own block of indices
@test_throws ArgumentError change_variance(Pv, Val((con2,cov2)))
@test_throws ArgumentError change_variance(Pv, Val((con2,cov2,con2)))

# A product encodes no index symmetry of its own, so a congruence of all its
# indices only gives their new dimension
pc = ProductTensorValue{CoVectorValue{cov2},CoVectorValue{cov2}}(1.0,2.0,3.0,4.0)
@test congruent_prod(pc, TensorValue{con2,con3}(1:6...)) isa TensorValue{con3,con3}
@test get_array(congruent_prod(pc, TensorValue{con2,con2}(1,2,3,4))) ≈
      get_array(congruent_prod(full_rank(pc), TensorValue{con2,con2}(1,2,3,4)))

# The variance of both factors is carried by the product size

pv = product_tensor(CoVectorValue(1.0,2.0), ContraVectorValue(1.0,2.0,3.0))
@test size(pv) === (cov2,con3)
@test tensor_type(pv) == (1,1,0)

# Components, conversions and comparison, against the dense `outer` product.
# Integer components keep the comparisons exact: expanding a dependent
# component of a factor distributes the product over a sum, which does not
# round like the dense product does.

for Ta in tags, Tb in tags
  A, B = make_concretetype(Ta{Int}), make_concretetype(Tb{Int})
  a = A(ntuple(identity, num_indep_components(A)))
  b = B(ntuple(i -> 10i, num_indep_components(B)))
  p = product_tensor(a,b)

  @test size(p) == (size(a)...,size(b)...)
  @test num_components(p) == num_components(a)*num_components(b)
  @test num_indep_components(p) == num_indep_components(a)*num_indep_components(b)
  # One stored component per pair of components stored by the factors, which
  # exceeds the independent ones for a factor that stores dependent components
  @test length(Tuple(p)) == length(Tuple(a))*length(Tuple(b))
  @test length(get_indep_components(p)) == num_indep_components(p)

  # The full tensor holds every independent component, so it reconstructs the
  # product, built from the canonical factor types that the product stores
  Ca, Cb = _minimal_tensortype_with_shape(typeof(a)), _minimal_tensortype_with_shape(typeof(b))
  @test ProductTensorValue{Ca,Cb}(get_array(p)) === p

  for cia in CartesianIndices(size(get_array(a))), cib in CartesianIndices(size(get_array(b)))
    @test p[Tuple(cia)...,Tuple(cib)...] == a[cia]*b[cib]
  end

  @test get_array(p) == get_array(outer(a,b))
  @test p ≈ outer(a,b)
  @test p == outer(a,b)
  @test outer(a,b) == p
  @test p == product_tensor(a,b)
  @test p != product_tensor(2a,b)

  # A product of tensors of a variance that cannot contract with itself has no
  # metric-free inner product, like the tensors themselves
  if all(i -> can_contract(i,i), size(p))
    @test inner(p,p) == inner(outer(a,b),outer(a,b))
  else
    @test_throws ArgumentError inner(p,p)
  end

  @test iszero(zero(p))
  @test p + zero(p) == p
  @test p - p == zero(p)

  # `component_basis` spans the product space and `indep_comp_getindex` reads
  # the coordinates in that basis
  basis = component_basis(typeof(p))
  @test length(basis) == num_indep_components(p)
  @test sum(indep_comp_getindex(p,i)*basis[i] for i in 1:num_indep_components(p); init=zero(p)) == p
end

# Random values of floating point components

for Ta in tags, Tb in tags
  A, B = make_concretetype(Ta{Float64}), make_concretetype(Tb{Float64})
  a, b = rand(A), rand(B)
  p = product_tensor(a,b)
  @test eltype(p) === Float64
  @test get_array(p) ≈ get_array(outer(a,b))
  @test p ≈ outer(a,b)
  @test 2p ≈ p + p
end

# A product of two symmetric tensors is a symmetric fourth order tensor, up to
# the component order of `SymFourthOrderTensorValue` (which stores the second
# factor's components contiguously)

a = rand(SymTensorValue{3,Float64})
b = rand(SymTensorValue{3,Float64})
p = product_tensor(a,b)
s4 = outer(a,b)
@test isa(s4,SymFourthOrderTensorValue{3})
@test get_array(p) ≈ get_array(s4)
@test p ≈ s4
@test Tuple(p) != Tuple(s4)

# Conversions from and to the full rank tensor of same size
#
# A full rank factor holds one component per entry of its full tensor, in the
# column major order of that tensor.

function _holds_full_tensor(V)
  s = Int.(size(V))
  num_indep_components(V) == prod(s) &&
    _indep_comp_cartesian_indices(V) == [Tuple(I) for I in vec(CartesianIndices(s))]
end

@test all(V -> is_full_rank(V) === _holds_full_tensor(V), tags)

# Each stored component, dependent ones included, is read from the full tensor
# at its `_stored_comp_cartesian_indices`, in the order they are stored

function _reads_stored_comps(V)
  VV = make_concretetype(_minimal_tensortype_with_shape(V){Int})
  u = VV(ntuple(i -> i, num_indep_components(VV)))
  Is = _stored_comp_cartesian_indices(VV)
  length(Is) == num_stored_components(VV) &&
    all(get_array(u)[I...] == Tuple(u)[k] for (k,I) in enumerate(Is))
end

@test all(_reads_stored_comps, tags)
@test _reads_stored_comps(SymFourthOrderTensorValue{2})
# a type that stores dependent components needs its own method, and gives more
# indices than it has independent components
@test length(_stored_comp_cartesian_indices(SymTracelessTensorValue{3})) == 6
@test_throws ErrorException _stored_comp_cartesian_indices(
  ProductTensorValue{VectorValue{2},SymTracelessTensorValue{2},Int})

# The slot map holds `±k` at every full tensor index that reads the stored
# component `k`, and `0` where the full tensor vanishes identically, so every
# component of a tensor is read from its full tensor through it. This is what the
# weighted contraction below rests on: a type whose full tensor held a
# combination of several stored components would break it silently.

function _reads_through_slot_map(V)
  VV = make_concretetype(_minimal_tensortype_with_shape(V){Float64})
  slots = get_array(stored_comp_indices(VV))
  all(1:20) do _
    u = VV(ntuple(i -> randn(), num_indep_components(VV)))
    full, stored = get_array(u), Tuple(u)
    all(CartesianIndices(slots)) do I
      k = slots[I]
      iszero(k) ? iszero(full[I]) : full[I] == sign(k)*stored[abs(k)]
    end
  end
end

@test all(_reads_through_slot_map, tags)
@test all(_reads_through_slot_map, (SymFourthOrderTensorValue{2}, SymFourthOrderTensorValue{3},
                                    SkewSymTensorValue{4}, SymNTensorValue{Tuple{3,3,3}}))

# A type materialising a component that no single index reads has no slot map

@test_throws ErrorException stored_comp_indices(
  ProductTensorValue{VectorValue{2},SymTracelessTensorValue{2},Int})

# A symmetric factor materialises no such component, so the product keeps a map

@test Tuple(stored_comp_indices(
  ProductTensorValue{VectorValue{2},SymTensorValue{2},Int})) == ntuple(identity,6)

# The weights count the indices reading each stored component, that is the
# Euclidean multiplicities of the stored layout. They do not depend on the
# variance, and a traceless tensor shares the symmetric layout's

@test _stored_comp_weights(SymTensorValue{3})             == [1,2,2,1,2,1]
@test _stored_comp_weights(SymTracelessTensorValue{3})    == [1,2,2,1,2,1]
@test _stored_comp_weights(SymTensorValue{con3})          == [1,2,2,1,2,1]
@test _stored_comp_weights(SkewSymTensorValue{3})         == [2,2,2]
@test _stored_comp_weights(SymNTensorValue{Tuple{3,3,3}}) == [1,3,3,3,6,3,1,3,3,1]
@test_throws ErrorException _stored_comp_weights(
  ProductTensorValue{VectorValue{2},SymTracelessTensorValue{2},Int})

# An alternating K-tensor stores one component per sorted combination of indices,
# read at each of its `K!` permutations, and none at all when `K > D`

@test _stored_comp_weights(ExteriorFormValue{1,cov3})   == fill(1, 3)
@test _stored_comp_weights(ExteriorFormValue{2,cov3})   == fill(2, 3)
@test _stored_comp_weights(ExteriorFormValue{3,Cov(4)}) == fill(6, 4)
@test isempty(_stored_comp_weights(ExteriorFormValue{4,cov3}))

# Two tensors sharing a slot map pair their components up with cancelling signs,
# so their Euclidean contraction is the sum over the stored components weighted
# by those multiplicities, and agrees with the dense contraction over the full
# tensors. Sharing the map, and not only the multiplicities, is what makes it
# hold for a symmetric and a traceless tensor against each other

s3, t3 = SymTensorValue{3}(1.,2,3,4,5,6), SymTracelessTensorValue{3}(1.,2,3,4,5)
@test flat_inner(s3, t3) == 63.0
@test flat_inner(s3, t3) == _dense_flat_inner(s3, t3)
@test flat_inner(t3, s3) == _dense_flat_inner(t3, s3)
@test flat_inner(t3, t3) == _dense_flat_inner(t3, t3)
@test flat_inner(s3, s3) == _dense_flat_inner(s3, s3)

@test flat_inner(SkewSymTensorValue{3}(1.,2,3), SkewSymTensorValue{3}(1.,2,3)) == 28.0
for u in (SymNTensorValue{Tuple{3,3,3}}(ntuple(Float64, 10)),
          SymFourthOrderTensorValue{3}(ntuple(Float64, 36)),
          ExteriorFormValue{2,cov3}(1.,2,3),
          ExteriorVectorValue{3,Con(4)}(1.,2,3,4))
  @test flat_inner(u, u) == _dense_flat_inner(u, u)
end

# A product contracts through the same weights, whether or not its factors carry
# variance: the multiplicities are Euclidean

for p in (ProductTensorValue{VectorValue{2},SymTensorValue{2}}(ntuple(Float64, 6)),
          ProductTensorValue{ContraVectorValue{con2},SymTensorValue{con2}}(ntuple(Float64, 6)),
          ProductTensorValue{SymTensorValue{con3},SymTensorValue{con3}}(ntuple(Float64, 36)))
  @test flat_inner(p, p) == _dense_flat_inner(p, p)
end

# The predicate is constant folded, so the conversions carry no run time check

_full_rank_val(::Type{V}) where V = Val(is_full_rank(V))
@test @inferred(_full_rank_val(TensorValue{cov2,con3})) === Val(true)
@test @inferred(_full_rank_val(CoVectorValue{cov2})) === Val(true)
@test @inferred(_full_rank_val(SymTensorValue{3})) === Val(false)

@test is_full_rank(VectorValue{3})
@test is_full_rank(TensorValue{2,3})
@test is_full_rank(ThirdOrderTensorValue{2,3,4})
@test is_full_rank(CoVectorValue{cov2})
@test is_full_rank(ContraVectorValue{con3})
@test !is_full_rank(SymTensorValue{2})
@test !is_full_rank(ExteriorFormValue{2,cov3})

"""
    _full_rank_tensortype(::Type{<:ProductTensorValue{A,B,T,S,N}})

The full rank tensor type holding the same components as the product type: a
`TensorValue` for an order two product, an `NTensorValue` beyond. Its size
parameter is the product's, so the variance of each dimension is preserved.

`A` and `B` must be full rank for both types to hold the same components, which
this function does not check.
"""
@generated function _full_rank_tensortype(
  ::Type{<:ProductTensorValue{A,B,T,S,N}}) where {A,B,T,S,N}

  s = _Size_to_size(S)
  L = prod(_rm_size_variance(s))
  V = N == 2 ? TensorValue{s...,T,L} :
      N  > 2 ? NTensorValue{S,T,N,L} :
      throw(ErrorException("A product of tensors of order at least one has order at least two, got order $N"))
  :( $V )
end

full_tags = filter(is_full_rank, tags)

for Ta in full_tags, Tb in full_tags
  A, B = make_concretetype(Ta{Int}), make_concretetype(Tb{Int})
  a = A(ntuple(identity, num_indep_components(A)))
  b = B(ntuple(i -> 10i, num_indep_components(B)))
  p = product_tensor(a,b)
  F = _full_rank_tensortype(typeof(p))
  t = convert(F,p)

  @test t isa F
  @test size(t) === size(p)
  # the data field is passed through untouched, in both directions, which is what
  # makes the conversion cost nothing
  @test t.data === p.data
  @test convert(typeof(p),t).data === p.data
  # and the full tensors agree, which is what a wrong component order would break
  @test get_array(t) == get_array(p)
  @test t == p

  @test convert(typeof(p),t) === p
  @test convert(ProductTensorValue{_minimal_tensortype_with_shape(A),_minimal_tensortype_with_shape(B),Int},t) === p
  @test convert(change_eltype(F,Float64),p) == t

  # the product is the common type of the pair: it holds the same components and
  # knows the two factors on top, so arithmetic with the full rank tensor lands on it
  @test promote_type(typeof(p), F) === typeof(p)
  @test promote_type(F, typeof(p)) === typeof(p)
  @test p + t === p + p
end

pv = product_tensor(CoVectorValue(1.0,2.0), ContraVectorValue(1.0,2.0,3.0))
Fv = _full_rank_tensortype(typeof(pv))
@test Fv === TensorValue{cov2,con3,Float64,6}
@test typeof(@inferred(convert(Fv,pv))) === Fv
@test @inferred(convert(typeof(pv),convert(Fv,pv))) === pv

pw = product_tensor(TensorValue{cov2,con3,Float64}(ntuple(Float64,6)), CoVectorValue(1.0,2.0))
@test _full_rank_tensortype(typeof(pw)) === ThirdOrderTensorValue{cov2,con3,cov2,Float64,12}

pe = product_tensor(zero(VectorValue{0,Float64}), VectorValue(1.0,2.0))
@test _full_rank_tensortype(typeof(pe)) === TensorValue{0,2,Float64,0}
@test Tuple(convert(TensorValue{0,2,Float64,0},pe)) === ()
@test convert(typeof(pe),convert(TensorValue{0,2,Float64,0},pe)) === pe

p3 = product_tensor(rand(ThirdOrderTensorValue{2,2,2,Float64}), VectorValue(1.0,2.0))
F3 = _full_rank_tensortype(typeof(p3))
@test F3 === NTensorValue{Tuple{2,2,2,2},Float64,4,16}
@test get_array(convert(F3,p3)) == get_array(p3)
@test convert(typeof(p3),convert(F3,p3)) === p3

# Each direction is one method taking any tensor type of the product's size, so
# a partially parametrized target is enough, and a product of other factors
# holding the same components converts like any other tensor

@test @inferred(convert(ProductTensorValue{CoVectorValue{cov2},ContraVectorValue{con3},Float64},convert(Fv,pv))) === pv
@test typeof(@inferred(convert(NTensorValue{Tuple{2,2,2,2},Float64},p3))) === F3

pp = convert(ProductTensorValue{ScalarValue,TensorValue{cov2,con3},Float64}, pv)
@test pp.data === pv.data
@test get_array(pp) == get_array(pv)
@test convert(typeof(pv),pp) === pv

# A `ScalarValue` factor adds no dimension, so a product holding one is a tensor
# of the size of its other factor, and converts with it. Two of them make a
# product of size (), the one component tensor that `ScalarValue` also is.

pvs = product_tensor(VectorValue(1.0,2.0), ScalarValue(3.0))
@test typeof(pvs) === ProductTensorValue{VectorValue{2},ScalarValue,Float64,Tuple{2},1,2}
@test @inferred(convert(VectorValue{2,Float64},pvs)) === VectorValue(3.0,6.0)
@test @inferred(convert(typeof(pvs),VectorValue(3.0,6.0))) === pvs

pss = product_tensor(ScalarValue(2.0), ScalarValue(3.0))
@test typeof(pss) === ProductTensorValue{ScalarValue,ScalarValue,Float64,Tuple{},0,1}
@test @inferred(convert(ScalarValue{Float64},pss)) === ScalarValue(6.0)
@test @inferred(convert(typeof(pss),ScalarValue(6.0))) === pss
@test convert(ProductTensorValue{ScalarValue,ScalarValue,Float64},ScalarValue(6.0)) === pss
@test promote_type(typeof(pss),ScalarValue{Float64}) === typeof(pss)
@test pss + ScalarValue(1.0) === convert(typeof(pss),ScalarValue(7.0))

# The other factor keeps its own components, whether it is full rank or not, and
# only its own type reads them: no other tensor of that size converts

psn = product_tensor(SymNTensorValue{Tuple{cov2,cov2},Float64}(1.0,2.0,3.0), ScalarValue(2.0))
SN = SymNTensorValue{Tuple{cov2,cov2},Float64,2,3}
@test typeof(psn) === ProductTensorValue{SymNTensorValue{Tuple{cov2,cov2}},ScalarValue,Float64,Tuple{cov2,cov2},2,3}
@test !is_full_rank(typeof(psn))
@test @inferred(convert(SN,psn)) === SN(2.0,4.0,6.0)
@test convert(SN,psn).data === psn.data
@test @inferred(convert(typeof(psn),SN(2.0,4.0,6.0))) === psn
@test convert(ProductTensorValue{SymNTensorValue{Tuple{cov2,cov2}},ScalarValue,Float64},SN(2.0,4.0,6.0)) === psn

@test_throws MethodError convert(TensorValue{cov2,cov2,Float64},psn)
@test_throws MethodError convert(typeof(psn),TensorValue{cov2,cov2,Float64}(ntuple(Float64,4)))

# The scalar factor is read in either position

pns = product_tensor(ScalarValue(2.0), SymNTensorValue{Tuple{cov2,cov2},Float64}(1.0,2.0,3.0))
@test @inferred(convert(SN,pns)) === SN(2.0,4.0,6.0)
@test @inferred(convert(typeof(pns),SN(2.0,4.0,6.0))) === pns

# A factor materialising its dependent components has the product materialise
# them too, so the data is passed through there as well

ptl = product_tensor(SymTracelessTensorValue(1.0,2.0,3.0,4.0,5.0), ScalarValue(2.0))
TL = SymTracelessTensorValue{3,Float64,6}
@test (num_indep_components(ptl), num_stored_components(typeof(ptl))) == (5,6)
@test @inferred(convert(TL,ptl)) === 2*SymTracelessTensorValue(1.0,2.0,3.0,4.0,5.0)
@test convert(TL,ptl).data === ptl.data
@test @inferred(convert(typeof(ptl),convert(TL,ptl))) === ptl

# Promotion stays on the full rank tensors, so these pairs have no common type
@test promote_type(typeof(psn),SN) !== typeof(psn)
@test promote_type(typeof(ptl),TL) !== typeof(ptl)

# Only a product of full rank factors converts: a factor holding fewer
# components than its full tensor has no full rank counterpart

ps = product_tensor(VectorValue(1,2), SymTensorValue(1,2,3))
@test_throws MethodError convert(NTensorValue{Tuple{2,2,2},Int}, ps)
@test_throws MethodError convert(
  ProductTensorValue{VectorValue{2},SymTensorValue{2},Int}, NTensorValue{Tuple{2,2,2},Int,3}(ntuple(identity,8)))

# `_full_rank_tensortype` reads the size only, it does not check the factors
@test _full_rank_tensortype(typeof(ps)) === NTensorValue{Tuple{2,2,2},Int,3,8}

# The full rank tensor must have the size of the product, variance included,
# which the two component counts alone do not tell apart. A tensor of another
# size has no conversion to the product at all

@test_throws MethodError convert(
  ProductTensorValue{VectorValue{2},VectorValue{3},Int}, TensorValue{3,2,Int}(ntuple(identity,6)))
@test_throws MethodError convert(
  ProductTensorValue{CoVectorValue{cov2},ContraVectorValue{con3},Int}, TensorValue{2,3,Int}(ntuple(identity,6)))
@test_throws MethodError convert(
  typeof(pv), TensorValue{con3,cov2,Float64,6}(ntuple(Float64,6)))

# Promotion with a full rank tensor keeps the factors and promotes the components

@test promote_type(change_eltype(typeof(pv),Int), Fv) === typeof(pv)
@test promote_type(Fv, change_eltype(typeof(pv),Int)) === typeof(pv)
@test promote_type(typeof(p3), F3) === typeof(p3)

tv = convert(Fv,pv)
@test pv + tv === pv + pv
@test tv - pv === zero(typeof(pv))
@test typeof(@inferred(pv + tv)) === typeof(pv)

# A full rank tensor of another size, variance included, is another tensor: no rule
# applies, the promotion falls back to the join of the two and the sum is rejected
@test promote_type(typeof(pv), TensorValue{con3,cov2,Float64,6}) !== typeof(pv)
@test promote_type(typeof(pv), TensorValue{2,3,Float64,6}) !== typeof(pv)
@test_throws ErrorException pv + TensorValue{con3,cov2,Float64,6}(ntuple(Float64,6))
@test_throws ErrorException pv + TensorValue{2,3,Float64,6}(ntuple(Float64,6))

# A tensor of the product's size but holding other components, because a factor
# holds fewer than its full tensor, has no common type with it either: the two
# would not promote to one another's components
@test promote_type(typeof(ps), NTensorValue{Tuple{2,2,2},Int,3,8}) !== typeof(ps)
@test promote_type(NTensorValue{Tuple{2,2,2},Int,3,8}, typeof(ps)) !== typeof(ps)
@test_throws MethodError ps + NTensorValue{Tuple{2,2,2},Int,3}(ntuple(identity,8))

# Two products of different factors have no common type, although they convert
# to one another: neither of the two factorisations holds the other
@test promote_type(typeof(pv), typeof(pp)) !== typeof(pv)
@test promote_type(typeof(pv), typeof(pp)) !== typeof(pp)
@test promote_type(typeof(pv), change_eltype(typeof(pv),Int)) === typeof(pv)

# VTK component names follow the storage order

p = product_tensor(VectorValue(1,2), SymTensorValue(1,2,3))
@test indep_components_names(typeof(p)) == ["XXX","YXX","XXY","YXY","XYY","YYY"]

# Empty factors

p = zero(ProductTensorValue{VectorValue{0},SymTensorValue{2},Float64,Tuple{0,2,2},3,0})
@test num_indep_components(p) == 0
@test Tuple(p) == ()
@test iszero(inner(p,p))
@test_throws ErrorException ProductTensorValue{VectorValue{0},SymTensorValue{2}}(())

# The factors give the component type that the empty data cannot
@test typeof(ProductTensorValue{VectorValue{0,Float64},SymTensorValue{2}}(())) === typeof(p)
@test typeof(product_tensor(zero(VectorValue{0,Float64}),zero(SymTensorValue{2,Int}))) === typeof(p)
@test typeof(ProductTensorValue{VectorValue{0},SymTensorValue{2},Float64}(())) === typeof(p)

# Order 0 factors

# A ScalarValue brings no index, so the product has the shape and the components
# of its other factor alone: this is what keeps the derivation indices of the
# partial derivatives of a scalar field in a product
p0 = product_tensor(VectorValue(1.0,2.0), ScalarValue(3.0))
@test typeof(p0) === ProductTensorValue{VectorValue{2},ScalarValue,Float64,Tuple{2},1,2}
@test size(p0) === (2,)
@test Tuple(p0) === (3.0,6.0)
@test get_array(p0) == get_array(VectorValue(3.0,6.0))
@test indep_components_names(typeof(p0)) == ["X","Y"]

# The canonical order 0 tag holds no component type, the data gives it
@test typeof(product_tensor(VectorValue(1,2), ScalarValue(3.0))) === typeof(p0)

# In first position it leaves the counts and the storage layout of the other
# factor untouched, both components counts included
ps0 = product_tensor(ScalarValue(2.0), SymTensorValue(1.0,2.0,3.0))
@test typeof(ps0) === ProductTensorValue{ScalarValue,SymTensorValue{2},Float64,Tuple{2,2},2,3}
@test (num_indep_components(ps0), num_stored_components(ps0)) == (3,3)
@test get_array(ps0) == get_array(2*SymTensorValue(1.0,2.0,3.0))

# including a factor that materialises its dependent components
pt0 = product_tensor(ScalarValue(2.0), SymTracelessTensorValue(1.0,2.0))
@test (num_indep_components(pt0), num_stored_components(pt0)) == (2,3)
@test get_array(pt0) == get_array(2*SymTracelessTensorValue(1.0,2.0))

# Wrong number of components

@test_throws AssertionError V(1.0,2.0)

# A product is not a valid factor, products do not nest

pp = product_tensor(VectorValue(1,2),VectorValue(3,4))
PP = typeof(pp)
@test_throws ErrorException _minimal_tensortype_with_shape(PP)
@test_throws ErrorException product_tensor(pp,VectorValue(1,2))
@test_throws ErrorException product_tensor(VectorValue(1,2),pp)
@test_throws ErrorException product_tensor(pp,pp)
@test_throws ErrorException ProductTensorValue{PP,VectorValue{2}}(ntuple(identity,8))
@test_throws ErrorException ProductTensorValue{PP,VectorValue{2},Int}(ntuple(identity,8))
@test_throws ErrorException ProductTensorValue{PP,VectorValue{2}}(())

# A product contraction consumes a whole factor of each operand: the second
# factor of the first operand against the first factor of the second one

pc_a = ProductTensorValue{SymTensorValue{2},ContraVectorValue{con2}}(ntuple(Float64, 6))
pc_b = ProductTensorValue{CoVectorValue{cov2},VectorValue{3}}(ntuple(Float64, 6))
pc = dot(pc_a, pc_b)
@test pc isa ProductTensorValue{SymTensorValue{2},VectorValue{3}}
@test get_array(pc) ≈ [ sum(get_array(pc_a)[i,j,l]*get_array(pc_b)[l,k] for l in 1:2)
                        for i in 1:2, j in 1:2, k in 1:3 ]

# A factor of order 2 leaves no single index for the contraction to consume

@test_throws ErrorException dot(
  ProductTensorValue{SymTensorValue{2},SymTensorValue{2}}(ntuple(Float64, 9)), pc_b)

# Point as ContraVectorValue overloads

pp_a   = ProductTensorValue{SymNTensorValue{Tuple{cov2,cov2}},VectorValue{2}}(ntuple(Float64, 6))
pp_flat = ProductTensorValue{VectorValue{2},VectorValue{3}}(ntuple(Float64, 6))
pp_ref  = [ sum(get_array(pp_a)[i,j,l]*get_array(pp_flat)[l,k] for l in 1:2)
            for i in 1:2, j in 1:2, k in 1:3 ]

@test dot(pp_a, pp_flat) isa ProductTensorValue{SymNTensorValue{Tuple{cov2,cov2}},VectorValue{3}}
@test get_array(dot(pp_a, pp_flat)) == pp_ref
@test get_array(dot(pp_a, pc_b)) == pp_ref
# Int - Con is wrong
@test_throws DimensionMismatch dot(
  pp_a, ProductTensorValue{ContraVectorValue{con2},VectorValue{3}}(ntuple(Float64, 6)))

# A partials tensor keeps its derivation indices through the operations acting on
# its second factor alone

∂u  = ProductTensorValue{SymNTensorValue{Tuple{cov2}},ScalarValue}(1.0,2.0)
∂∂u = ProductTensorValue{SymNTensorValue{Tuple{cov2,cov2}},ScalarValue}(1.0,2.0,3.0)
∂v  = ProductTensorValue{SymNTensorValue{Tuple{cov2}},ContraVectorValue{con2}}(ntuple(Float64,4))
ω   = CoVectorValue{cov2}(1.0,2.0)
x   = VectorValue(1.0,2.0)

@test outer(∂u, ω) isa ProductTensorValue{SymNTensorValue{Tuple{cov2}},CoVectorValue{cov2}}
@test get_array(outer(∂∂u, ω)) ≈ get_array(contracted_product(Val(0), ∂∂u, ω))
@test dot(∂v, ω) isa ProductTensorValue{SymNTensorValue{Tuple{cov2}},ScalarValue}
@test get_array(dot(∂v, ω)) ≈ get_array(contracted_product(Val(1), ∂v, ω))

# The second factor of the right operand is kept factored when a whole factor of
# each is consumed, and flattened when the contraction has to reach inside one
∂a = ProductTensorValue{SymNTensorValue{Tuple{cov2}},SymTensorValue{con2}}(ntuple(Float64,6))
@test dot(∂v, ∂a) isa ProductTensorValue{SymNTensorValue{Tuple{cov2}},SymTensorValue{con2}}
@test dot(∂v, product_tensor(SymNTensorValue{Tuple{cov2,cov2}}(1.0,2.0,3.0), ω)) isa
      ProductTensorValue{SymNTensorValue{Tuple{cov2}},TensorValue{cov2,cov2}}

# The contraction reaches the last derivation index when the second factor has no
# index of its own, so there the derivation indices cannot be kept apart
@test dot(∂∂u, x) isa CoVectorValue{cov2}
@test get_array(dot(∂∂u, x)) ≈ get_array(∂∂u)*get_array(x)
@test dot(∂∂u, ContraVectorValue{con2}(1.0,2.0)) === dot(∂∂u, x)

# A factor with no independent component has no slice to act on, so the dense
# product is the only possible result
∂w = product_tensor(zero(SymNTensorValue{Tuple{cov2},Float64}),
                    zero(ExteriorFormValue{3,cov2,Float64}))
@test outer(∂w, x) isa NTensorValue{Tuple{cov2,cov2,cov2,cov2,con2}}
@test iszero(outer(∂w, x))

# A point coefficient contributes a contravariant index, so that the jacobian it
# builds has opposite variances
@test outer(∂u, x) isa ProductTensorValue{SymNTensorValue{Tuple{cov2}},ContraVectorValue{con2}}
@test get_array(outer(∂u, x)) ≈ get_array(∂u)*get_array(x)'
@test is_opposit_variance(size(outer(∂u, x))...)
@test outer(∂u, x) === outer(∂u, ContraVectorValue{con2}(1.0,2.0))

end # module
