module bm_tensor_values_microbench

using PkgBenchmark, BenchmarkTools
using Gridap, Gridap.TensorValues

SUITE = BenchmarkGroup()

# Tensor operations, grouped by operation and dimension. The labels give the order of the operands.
# The result type of each operation is checked in test/TensorAlgebraTests/TensorTypes.jl.
# The operations which are not supported yet are skipped, and they are benchmarked once they are implemented.

# The operands are passed through a `Ref` to avoid constant propagation
function add_tensor_benchmark!(group, label, f, args...)
  try
    f(args...)
  catch
    return group  # not supported yet
  end
  if length(args) == 1
    group[label] = @benchmarkable $f($(Ref(args[1]))[])
  else
    group[label] = @benchmarkable $f($(Ref(args[1]))[], $(Ref(args[2]))[])
  end
  group
end

function add_tensor_operations_benchmarks!(SUITE, D)
  v = VectorValue{D,Float64}(ntuple(i -> Float64(i), D))
  u = VectorValue{D,Float64}(ntuple(i -> Float64(D + i), D))
  A = TensorValue{D,D,Float64}(ntuple(i -> 1.0 + i / 10, D^2))
  B = TensorValue{D,D,Float64}(ntuple(i -> 2.0 - i / 10, D^2))
  W = ThirdOrderTensorValue{D,D,D,Float64}(ntuple(i -> 1.0 + i / 100, D^3))
  H = HighOrderTensorValue{Tuple{D,D,D,D},Float64}(ntuple(i -> 1.0 + i / 1000, D^4))
  dim = "$(D)D"

  group = SUITE["TensorValues"]["outer"][dim]
  for (label, a, b) in (("1⊗1", v, u), ("1⊗2", v, A), ("2⊗1", A, v), ("2⊗2", A, B),
                        ("1⊗3", v, W), ("3⊗1", W, v))
    add_tensor_benchmark!(group, label, ⊗, a, b)
  end

  group = SUITE["TensorValues"]["inner"][dim]
  for (label, a, b) in (("1⊙1", v, u), ("2⊙2", A, B), ("3⊙3", W, W), ("4⊙4", H, H))
    add_tensor_benchmark!(group, label, ⊙, a, b)
  end

  group = SUITE["TensorValues"]["dot"][dim]
  for (label, a, b) in (("1⋅1", v, u), ("1⋅2", v, A), ("2⋅1", A, v), ("2⋅2", A, B),
                        ("1⋅3", v, W), ("3⋅1", W, v), ("2⋅3", A, W), ("3⋅2", W, A), ("3⋅3", W, W),
                        ("2⋅4", A, H), ("4⋅2", H, A))
    add_tensor_benchmark!(group, label, ⋅, a, b)
  end

  group = SUITE["TensorValues"]["double dot"][dim]
  for (label, a, b) in (("2⋅²2", A, B), ("4⋅²2", H, A), ("2⋅²4", A, H), ("4⋅²4", H, H), ("3⋅²2", W, A), ("3⋅²4", W, H))
    add_tensor_benchmark!(group, label, ⋅², a, b)
  end

  if D == 3
    group = SUITE["TensorValues"]["cross"][dim]
    for (label, a, b) in (("1×1", v, u),)# ("2×2", A, B), ("2×4", A, H), ("4×2", H, A), ("3×2", W, A))
      add_tensor_benchmark!(group, label, ×, a, b)
    end
  end

  group = SUITE["TensorValues"]["arithmetic"][dim]
  for (name, X) in (("1", v), ("2", A), ("3", W), ("4", H))
    add_tensor_benchmark!(group, "$name+$name", +, X, X)
    add_tensor_benchmark!(group, "$name-$name", -, X, X)
    add_tensor_benchmark!(group, "+$name", +, X)
    add_tensor_benchmark!(group, "-$name", -, X)
    add_tensor_benchmark!(group, "α*$name", *, 2.0, X)
    add_tensor_benchmark!(group, "$name/α", /, X, 2.0)
    add_tensor_benchmark!(group, "$name'", transpose, X)
  end

end

add_tensor_operations_benchmarks!(SUITE, 3)

end
