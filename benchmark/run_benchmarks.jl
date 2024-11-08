using Pkg

Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
Pkg.instantiate()

using Gridap
using PkgBenchmark

config_kwargs = (;
  juliacmd = `julia -O3`,
)

if haskey(ENV,"BM_TARGET") # Provided by CI workflow
  target = BenchmarkConfig(;config_kwargs..., id = ENV["BM_TARGET"])
else # Default to the current commit
  target = BenchmarkConfig(;config_kwargs...)
end

if haskey(ENV,"BM_BASE") # Provided by CI workflow
  base = BenchmarkConfig(;config_kwargs..., id = ENV["BM_BASE"])
else # Default to master
  base = BenchmarkConfig(;config_kwargs..., id = "master")
end

results = judge(Gridap, target, base)
outfile = normpath(@__DIR__,"benchmark_results.md")
export_markdown(outfile,results)
