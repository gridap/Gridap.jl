module ProductTensorValueTests

using Gridap.TensorValues
using Gridap.TensorValues: _minimal_tensortype_with_shape, _factor_eltype, stored_comp_indices
using Gridap.TensorValues: _prod_counts, _prod_eltype
using Gridap.TensorValues: get_indep_components, num_stored_components
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
  HighOrderTensorValue{Tuple{2,2,2}},
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
P = ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2},3,6}
@test typeof(ProductTensorValue(a,b)) === P

# A component type given in a factor is not kept in the factor tags, it gives
# the component type of the product instead

@test typeof(@inferred(ProductTensorValue{VectorValue{2,Float64},SymTensorValue{2}}(ntuple(Int,6)))) === P
@test typeof(ProductTensorValue{VectorValue{2},SymTensorValue{2,Float64,3}}(ntuple(Int,6))) === P
@test typeof(ProductTensorValue{VectorValue{2,Int},SymTensorValue{2,Float64,3}}(ntuple(Int,6))) === P
@test typeof(ProductTensorValue{VectorValue{2},SymTensorValue{2}}(ntuple(Float64,6))) === P
@test typeof(ProductTensorValue{VectorValue{2},SymTensorValue{2}}(ntuple(Int,6))) === change_eltype(P,Int)

# An explicitly given component type prevails over that of the factors, as long
# as the latter promotes to it

@test typeof(ProductTensorValue{VectorValue{2,Int},SymTensorValue{2,Int,3},Float64,Tuple{2,2,2},3,6}(ntuple(Int,6))) === P
@test_throws AssertionError ProductTensorValue{VectorValue{2,Float64},SymTensorValue{2},Int,Tuple{2,2,2},3,6}(ntuple(Int,6))

# `S`, `N` and `L` follow from the factor types, so the component type is the
# only other parameter to give

@test typeof(@inferred(ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64}(ntuple(Int,6)))) === P
@test typeof(ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64}(ntuple(Int,6)...)) === P
@test_throws AssertionError ProductTensorValue{VectorValue{2,Float64},SymTensorValue{2},Int}(ntuple(Int,6))

# `S`, `N` and `L` may still be spelled out, in part or in full, and are then
# checked against the ones the factors give

@test typeof(@inferred(ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2}}(ntuple(Int,6)))) === P
@test typeof(ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2}}(ntuple(Int,6)...)) === P
@test typeof(@inferred(ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2},3}(ntuple(Int,6)))) === P
@test typeof(ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2},3}(ntuple(Int,6)...)) === P

@test_throws AssertionError ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,3}}(ntuple(Int,6))
@test_throws AssertionError ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2},2}(ntuple(Int,6))

# All the spellings build the same value

@test allequal((
  ProductTensorValue{VectorValue{2},SymTensorValue{2}}(ntuple(Float64,6)),
  ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64}(ntuple(Int,6)),
  ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64}(ntuple(Int,6)...),
  ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2}}(ntuple(Int,6)),
  ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2}}(ntuple(Int,6)...),
  ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2},3}(ntuple(Int,6)),
  ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2},3}(ntuple(Int,6)...),
  ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2},3,6}(ntuple(Int,6)),
  ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2},3,6}(ntuple(Int,6)...)))

# Storage layout: the independent components of the first factor vary fastest,
# so a product of two tensors without dependent components has the components
# of their `outer` product, in the same order.

p = ProductTensorValue(VectorValue(1,2), VectorValue(10,20,30))
@test Tuple(p) == (10,20,20,40,30,60)
@test Tuple(p) == Tuple(VectorValue(1,2) ⊗ VectorValue(10,20,30))
@test size(p) == (2,3)
@test num_indep_components(p) == 6

p = ProductTensorValue(VectorValue(1,2), SymTensorValue(10,20,30))
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

pq = ProductTensorValue(q,v)
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

p0 = ProductTensorValue(zero(SymTracelessTensorValue{1,Float64}), v)
@test num_indep_components(p0) == 0
@test length(Tuple(p0)) == 2
@test eltype(p0) === Float64
@test iszero(p0)

# Type parameters and inference

a, b = VectorValue(1.0,2.0), SymTensorValue(1.0,2.0,3.0)
p = @inferred ProductTensorValue(a,b)
@test typeof(p) === ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2},3,6}
@test @inferred(zero(typeof(p))) == ProductTensorValue(zero(a),b)
@test typeof(@inferred(p+p)) === typeof(p)
@test typeof(@inferred(2.0*p)) === typeof(p)
@test typeof(@inferred(2*p)) === typeof(p)
@test @inferred(p[1,2,2]) == a[1]*b[2,2]
@test typeof(@inferred(ProductTensorValue{VectorValue{2},SymTensorValue{2},Float64,Tuple{2,2,2},3,6}(ntuple(Float64,6)))) === typeof(p)

@test change_eltype(p,Int) === ProductTensorValue{VectorValue{2},SymTensorValue{2},Int,Tuple{2,2,2},3,6}
@test promote_type(typeof(p),change_eltype(p,Int)) === typeof(p)
@test convert(change_eltype(p,Float32),p) == p

# The variance of both factors is carried by the product size

pv = ProductTensorValue(CoVectorValue(1.0,2.0), ContraVectorValue(1.0,2.0,3.0))
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
  p = ProductTensorValue(a,b)

  @test size(p) == (size(a)...,size(b)...)
  @test num_components(p) == num_components(a)*num_components(b)
  @test num_indep_components(p) == num_indep_components(a)*num_indep_components(b)
  # One stored component per pair of components stored by the factors, which
  # exceeds the independent ones for a factor that stores dependent components
  @test length(Tuple(p)) == length(Tuple(a))*length(Tuple(b))
  @test length(get_indep_components(p)) == num_indep_components(p)

  for cia in CartesianIndices(size(get_array(a))), cib in CartesianIndices(size(get_array(b)))
    @test p[Tuple(cia)...,Tuple(cib)...] == a[cia]*b[cib]
  end

  @test get_array(p) == get_array(outer(a,b))
  @test p ≈ outer(a,b)
  @test p == outer(a,b)
  @test outer(a,b) == p
  @test p == ProductTensorValue(a,b)
  @test p != ProductTensorValue(2a,b)

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
  p = ProductTensorValue(a,b)
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
p = ProductTensorValue(a,b)
s4 = outer(a,b)
@test isa(s4,SymFourthOrderTensorValue{3})
@test get_array(p) ≈ get_array(s4)
@test p ≈ s4
@test Tuple(p) != Tuple(s4)

# VTK component names follow the storage order

p = ProductTensorValue(VectorValue(1,2), SymTensorValue(1,2,3))
@test indep_components_names(typeof(p)) == ["XXX","YXX","XXY","YXY","XYY","YYY"]

# Empty factors

p = zero(ProductTensorValue{VectorValue{0},SymTensorValue{2},Float64,Tuple{0,2,2},3,0})
@test num_indep_components(p) == 0
@test Tuple(p) == ()
@test iszero(inner(p,p))
@test_throws ErrorException ProductTensorValue{VectorValue{0},SymTensorValue{2}}(())

# The factors give the component type that the empty data cannot
@test typeof(ProductTensorValue{VectorValue{0,Float64},SymTensorValue{2}}(())) === typeof(p)
@test typeof(ProductTensorValue(zero(VectorValue{0,Float64}),zero(SymTensorValue{2,Int}))) === typeof(p)
@test typeof(ProductTensorValue{VectorValue{0},SymTensorValue{2},Float64}(())) === typeof(p)

# Wrong number of components

@test_throws AssertionError ProductTensorValue{VectorValue{2},SymTensorValue{2}}(1.0,2.0)

# A product is not a valid factor, products do not nest

pp = ProductTensorValue(VectorValue(1,2),VectorValue(3,4))
PP = typeof(pp)
@test_throws ErrorException _minimal_tensortype_with_shape(PP)
@test_throws ErrorException ProductTensorValue(pp,VectorValue(1,2))
@test_throws ErrorException ProductTensorValue(VectorValue(1,2),pp)
@test_throws ErrorException ProductTensorValue(pp,pp)
@test_throws ErrorException ProductTensorValue{PP,VectorValue{2}}(ntuple(identity,8))
@test_throws ErrorException ProductTensorValue{PP,VectorValue{2},Int}(ntuple(identity,8))
@test_throws ErrorException ProductTensorValue{PP,VectorValue{2}}(())

end # module
