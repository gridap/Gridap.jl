module DifferentialGeometryOperatorsTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.TensorValues: get_indep_components
using Gridap.Fields

using StaticArrays
using ForwardDiff

const cov2 = Cov(2)

_flat(a::Number)     = SVector(a)
_flat(a::MultiValue) = SVector(Tuple(get_array(a)))
_fullsize(a::Number)     = ()
_fullsize(a::MultiValue) = size(get_array(a))

function ref_partials(f, x0::SVector{D}, idxs::Tuple) where D
  isempty(idxs) && return _flat(f(Point(x0)))
  T = eltype(x0)
  e = SVector{D}(ntuple(k -> k == first(idxs) ? one(T) : zero(T), D))
  ForwardDiff.derivative(t -> ref_partials(f, x0 + t*e, Base.tail(idxs)), zero(T))
end

function test_partial_derivatives(f, x::Point{D}, Nmax) where D
  fx = f(x)
  sz = _fullsize(fx)
  for N in 1:Nmax
    pd = partial_derivatives(f,x,Val(N))
    @test typeof(pd) === derivatives_type(typeof(fx),x,Val(N))
    @test size(pd) === (ntuple(_->Cov(D),N)..., size(fx)...)
    @test (@inferred partial_derivatives(f,x,Val(N))) === pd
    for I in CartesianIndices(ntuple(_->D,N))
      r = ref_partials(f,get_array(x),Tuple(I))
      for (n,J) in enumerate(CartesianIndices(sz))
        @test pd[Tuple(I)...,Tuple(J)...] ≈ r[n] atol=1e-9
      end
    end
  end
end

x2 = Point(0.3,-0.7)
x3 = Point(0.3,-0.7,0.2)

for a in (1.0, 2, SymTensorValue{2}(1.,2.,3.), SymTracelessTensorValue{2}(1.,2.),
          TensorValue{2,2}(1.,2.,3.,4.), ExteriorFormValue{2,3}(1.,2.,3.),
          ThirdOrderTensorValue{2,2,2}(1.,2.,3.,4.,5.,6.,7.,8.))
  @test length(get_indep_components(a)) == num_indep_components(typeof(a))
  @test all(l -> indep_comp_getindex(a,l) === get_indep_components(a)[l],
            1:num_indep_components(typeof(a)))
end

@test get_indep_components(1.0) === (1.0,)

# scalar valued, cross checked against ForwardDiff's own gradient and hessian
fs(x) = sin(x[1])*exp(x[2]) + x[1]^3*x[2]

test_partial_derivatives(fs, x2, 2)
@test get_array(partial_derivatives(fs,x2,Val(1))) ≈ ForwardDiff.gradient(fs∘Point,get_array(x2))
@test get_array(partial_derivatives(fs,x2,Val(2))) ≈ ForwardDiff.hessian(fs∘Point,get_array(x2))

fi(x) = 3

@test partial_derivatives(fs,x2,Val(0)) === fs(x2)
@test partial_derivatives(fi,x2,Val(0)) === 3.0
@test typeof(partial_derivatives(fi,x2,Val(0))) === derivatives_type(Int,x2,Val(0))

# A vector valued function has the index layout of `gradient`, ∂u[i,j] = ∂u_j/∂xⁱ
fv(x) = VectorValue(sin(x[1])*x[2], exp(x[2])+x[1]^3)

@test get_array(partial_derivatives(fv,x2,Val(1))) ≈ get_array(gradient(fv,x2))
test_partial_derivatives(fv, x2, 2)

# Every tensor family, including those whose stored components are dependent
w(l,x) = sin((l+1)*x[1]) + l*x[2]^2 + x[1]*x[2]/(l+1)
mk(V,x) = V(ntuple(l -> w(l,x), Val(num_indep_components(V))))

for V in (VectorValue{2}, CoVectorValue{cov2}, TensorValue{2,2},
          SymTensorValue{2}, SymTracelessTensorValue{2},
          SkewSymTensorValue{2}, ExteriorFormValue{1,2},
          ExteriorFormValue{2,2}, ContraVectorValue{Con(2)})
  test_partial_derivatives(x -> mk(V,x), x2, 2)
end

for V in (VectorValue{3}, SymTensorValue{3}, ExteriorFormValue{2,3},
          ThirdOrderTensorValue{3,3,3})
  test_partial_derivatives(x -> mk(V,x), x3, 2)
end

# A value type without independent component has constant, zero components
f0(x) = zero(ExteriorFormValue{3,2,Float64,0})

@test iszero(partial_derivatives(f0,x2,Val(1)))
@test typeof(partial_derivatives(f0,x2,Val(1))) ===
  derivatives_type(ExteriorFormValue{3,2,Float64,0},x2,Val(1))

@test partial_derivatives(fv,Val(2))(x2) === partial_derivatives(fv,x2,Val(2))
@test partial_derivatives(fv)(x2) === partial_derivatives(fv,x2,Val(1))

# The operators consuming the partials are declared, not implemented
pdv = partial_derivatives(fv,x2,Val(1))

end # module
