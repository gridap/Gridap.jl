module ProductTensorValueTests

using Gridap.TensorValues
using Gridap.TensorValues: _minimal_tensortype_with_shape, _factor_eltype, stored_comp_indices
using Gridap.TensorValues: _prod_counts, _prod_eltype
using Gridap.TensorValues: get_indep_components, num_stored_components
using Gridap.TensorValues: _indep_comp_cartesian_indices, is_full_rank
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
  # inner product, like the tensors themselves
  if all(i -> can_contract(i,i), size(p))
    @test inner(p,p) == inner(outer(a,b),outer(a,b))
  else
    @test_throws DimensionMismatch inner(p,p)
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

  L = prod(Int.(S.parameters))
  V = N == 2 ? TensorValue{S.parameters...,T,L} :
      N  > 2 ? NTensorValue{S,T,N,L}            :
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

# Only a product of full rank factors converts: a factor holding fewer
# components than its full tensor has no full rank counterpart

ps = product_tensor(VectorValue(1,2), SymTensorValue(1,2,3))
@test_throws ErrorException convert(NTensorValue{Tuple{2,2,2},Int}, ps)
@test_throws ErrorException convert(
  ProductTensorValue{VectorValue{2},SymTensorValue{2},Int}, NTensorValue{Tuple{2,2,2},Int,3}(ntuple(identity,8)))

# `_full_rank_tensortype` reads the size only, it does not check the factors
@test _full_rank_tensortype(typeof(ps)) === NTensorValue{Tuple{2,2,2},Int,3,8}

# The full rank tensor must have the size of the product, variance included,
# which the two component counts alone do not tell apart

@test_throws AssertionError convert(
  ProductTensorValue{VectorValue{2},VectorValue{3},Int}, TensorValue{3,2,Int}(ntuple(identity,6)))
@test_throws AssertionError convert(
  ProductTensorValue{CoVectorValue{cov2},ContraVectorValue{con3},Int}, TensorValue{2,3,Int}(ntuple(identity,6)))

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

# The rule claims the product for any matching size; whether the factors are full
# rank is reported by the conversion that promotion then applies
@test promote_type(typeof(ps), NTensorValue{Tuple{2,2,2},Int,3,8}) === typeof(ps)
@test_throws ErrorException ps + NTensorValue{Tuple{2,2,2},Int,3}(ntuple(identity,8))

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

end # module
