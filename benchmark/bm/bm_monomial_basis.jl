module bm_monomial_basis

using PkgBenchmark, BenchmarkTools
using Gridap
using Gridap.Polynomials
using Gridap.TensorValues
using StaticArrays

################################################
# src/Polynomials/MonomialBasis.jl: _set_value_!
################################################

gradient_type = Gridap.Fields.gradient_type

_set_value! = Gridap.Polynomials._set_value!

function set_value_driver(f,T,D,x,n)
  k = 1
  s = one(T)
  for i in 1:n
    k = f(x,s,k)
  end
end

function set_value_benchmarkable(D, T, V, n)
  C = num_indep_components(V)
  x = zeros(V,n*C)
  return @benchmarkable set_value_driver($_set_value!,$T,$D,$x,$n)
end

##################################################
# src/Polynomials/ModalC0Bases.jl: _set_value_mc0!
##################################################

_set_value_mc0! = Gridap.Polynomials._set_value_mc0!

function set_value_mc0_driver(f,T,D,x,n)
  k = 1
  s = one(T)
  for i in 1:n
    k = f(x,s,k,2)
  end
end

function set_value_mc0_benchmarkable(D, T, V, n)
  C = num_indep_components(V)
  x = zeros(V,2*n*C)
  return @benchmarkable set_value_mc0_driver($_set_value_mc0!,$T,$D,$x,$n)
end

###################################################
# src/Polynomials/MonomialBasis.jl: _set_gradient!
###################################################

 _set_gradient! = Gridap.Polynomials. _set_gradient!

function set_gradient_driver(f,T,D,V,x,n)
  k = 1
  s = VectorValue{D,T}(ntuple(_->one(T),D))
  for i in 1:n
    k = f(x,s,k,V)
  end
end

function set_gradient_benchmarkable(D, T, V, n)
  C = num_indep_components(V)
  G = gradient_type(V, zero(Point{D,T}))
  x = zeros(G,n*C);
  return @benchmarkable set_gradient_driver($_set_gradient!,$T,$D,$V,$x,$n)
end

#####################################################
# src/Polynomials/ModalC0Bases.jl: _set_gradient_mc0!
#####################################################

 _set_gradient_mc0! = Gridap.Polynomials. _set_gradient_mc0!

function set_gradient_mc0_driver(f,T,D,V,x,n)
  k = 1
  s = VectorValue{D,T}(ntuple(_->one(T),D))
  for i in 1:n
    k = f(x,s,k,1,V)
  end
end

function set_gradient_mc0_benchmarkable(D, T, V, n)
  C = num_indep_components(V)
  G = gradient_type(V, zero(Point{D,T}))
  x = zeros(G,n*C);
  return @benchmarkable set_gradient_mc0_driver($_set_gradient_mc0!,$T,$D,$V,$x,$n)
end

#################################################
# src/Polynomials/MonomialBasis.jl: _evaluate_1d!
#################################################

_evaluate_1d! = Gridap.Polynomials._evaluate_1d!

function evaluate_1d_driver(f,order,D,v,x_vec)
  for x in x_vec
    f(v,x,order,D)
  end
end

function evaluate_1d_benchmarkable(D, T, V, n)
  n = Integer(n/50)
  order = num_indep_components(V)
  v = zeros(D,order+1);
  x = rand(MVector{n,T})
  return @benchmarkable evaluate_1d_driver($_evaluate_1d!,$order,$D,$v,$x)
end

################################################
# src/Polynomials/MonomialBasis.jl:_gradient_1d!
################################################

_gradient_1d! = Gridap.Polynomials._gradient_1d!

function gradient_1d_driver(f,order,D,v,x_vec)
  for x in x_vec
    f(v,x,order,D)
  end
end

function gradient_1d_benchmarkable(D, T, V, n)
  n = Integer(n/10)
  order = num_indep_components(V)
  v = zeros(D,order+1);
  x = rand(MVector{n,T})
  return @benchmarkable gradient_1d_driver($_gradient_1d!,$order,$D,$v,$x)
end

################################################
# src/Polynomials/MonomialBasis.jl:_hessian_1d!
################################################

_hessian_1d! = Gridap.Polynomials._hessian_1d!

function hessian_1d_driver(f,order,D,v,x_vec)
  for x in x_vec
    f(v,x,order,D)
  end
end

function hessian_1d_benchmarkable(D, T, V, n)
  n = Integer(n/10)
  order = num_indep_components(V)
  v = zeros(D,order+1);
  x = rand(MVector{n,T})
  return @benchmarkable hessian_1d_driver($_hessian_1d!,$order,$D,$v,$x)
end

#####################
# benchmarkable suite
#####################

const SUITE = BenchmarkGroup()

const benchmarkables = (
  set_value_benchmarkable,
  set_value_mc0_benchmarkable,
  set_gradient_benchmarkable,
  set_gradient_mc0_benchmarkable,
  evaluate_1d_benchmarkable,
  gradient_1d_benchmarkable,
  hessian_1d_benchmarkable
)

const dims=(1, 2, 3, 5, 8)
const n = 3000
const T = Float64

for benchable in benchmarkables
  for D in dims
    TV = [
      VectorValue{D,T},
      TensorValue{D,D,T,D*D},
      SymTensorValue{D,T,Integer(D*(D+1)/2)},
      SymTracelessTensorValue{D,T,Integer(D*(D+1)/2)}
    ]

    for V in TV
      if V == SymTracelessTensorValue{1,T,1} continue end # no dofs
      name = "monomial_basis_$(D)D_$(V)_$(benchable)"
      SUITE[name] = benchable(D, T, V, n)
    end
  end
end

end # module
