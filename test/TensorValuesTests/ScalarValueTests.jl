module ScalarValueTests

using Gridap.TensorValues
using Gridap.TensorValues: num_stored_components, get_indep_components
using Gridap.Arrays: get_array
using Test
using StaticArrays

# Construction

a = ScalarValue(2.0)
@test a isa ScalarValue{Float64}
@test a isa MultiValue{Tuple{},Float64,0,1}
@test a isa Number

# Every MultiValue also takes its independent components as a tuple, and an
# order 0 tensor is shaped like a zero dimensional array
@test a === ScalarValue{Float64}(2.0)
@test a === ScalarValue((2.0,))
@test a === ScalarValue{Float64}((2.0,))
@test a === ScalarValue(fill(2.0))
@test a === ScalarValue(StaticArrays.Scalar(2.0))
@test a === ScalarValue{Float64}(StaticArrays.Scalar(2.0))

# The component converts to T
@test ScalarValue{Float64}(2) === a

# An order 0 tensor holds a single component
@test_throws ErrorException ScalarValue((1.0,2.0))
@test_throws ErrorException ScalarValue{Float64}((1.0,2.0))

# A component is never a tensor in disguise
@test_throws ErrorException ScalarValue(VectorValue(1.0,2.0))
@test_throws ErrorException ScalarValue{Vector{Float64}}(([1.0],))

# Only a zero dimensional array has the shape of an order 0 tensor
@test_throws MethodError ScalarValue([1.0])

# Shape and component counts

@test size(a) == ()
@test length(a) == 1
@test num_components(a) == 1
@test num_indep_components(a) == 1
@test num_stored_components(a) == 1
@test eltype(a) === Float64
@test change_eltype(ScalarValue{Float64},Int) === ScalarValue{Int}

# Indexing gives the component, as for a 0-form

@test a[] === 2.0
@test a[1] === 2.0
@test a[CartesianIndex()] === 2.0
@test indep_comp_getindex(a,1) === 2.0
@test get_indep_components(a) === (2.0,)
@test Tuple(a) === (2.0,)
@test_throws AssertionError a[2]

# Conversions with the equivalent zero dimensional arrays

@test get_array(a) === StaticArrays.Scalar(2.0)
@test convert(SArray{Tuple{},Float64,0,1},a) === StaticArrays.Scalar(2.0)
@test convert(MArray{Tuple{},Float64,0,1},a) == StaticArrays.Scalar(2.0)
@test MultiValue(StaticArrays.Scalar(2.0)) === a
@test convert(ScalarValue{Float64},a) === a
@test convert(ScalarValue{Int},ScalarValue(2)) === ScalarValue(2)

# A single component tensor converts from and to its component, through the
# generic MultiValue conversions
@test convert(Float64,a) === 2.0
@test convert(ScalarValue{Float64},2.0) === a

# The additive neutral is a tensor, the multiplicative one is the bare
# component: it neutralises the only multiplication there is, by a scalar

@test zero(a) === ScalarValue(0.0)
@test zero(ScalarValue{Float64}) === ScalarValue(0.0)
@test one(a) === 1.0
@test a*one(a) === a
@test iszero(zero(a))
@test !iszero(a)

# Promotion

@test promote_type(ScalarValue{Int},ScalarValue{Float64}) === ScalarValue{Float64}
@test promote(ScalarValue(1),ScalarValue(2.0)) === (ScalarValue(1.0),ScalarValue(2.0))

# The component basis of a one component tensor is its unit

@test component_basis(ScalarValue{Float64}) == [ScalarValue(1.0)]
@test indep_components_names(ScalarValue{Float64}) == [""]

# The MultiValue arithmetic carries through, since it goes through the
# independent components

@test a + a === ScalarValue(4.0)
@test a - a === ScalarValue(0.0)
@test -a === ScalarValue(-2.0)
@test 2*a === ScalarValue(4.0)
@test a*2 === ScalarValue(4.0)
@test a/2 === ScalarValue(1.0)
@test a + 1.0 === ScalarValue(3.0)
@test a ≈ ScalarValue(2.0)
@test a == ScalarValue(2.0)

# Printing shows the component list, like the other tensors

@test sprint(show,a) == "(2.0,)"
@test sprint(show,MIME"text/plain"(),a) == "ScalarValue{Float64}(2.0,)"

end # module
