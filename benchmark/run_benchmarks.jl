using Pkg

Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
Pkg.instantiate()

using Gridap
using PkgBenchmark
using DrWatson

results = judge(
  Gridap, 
  BenchmarkConfig(juliacmd = `julia -O3`), # target -> current branch
  BenchmarkConfig(juliacmd = `julia -O3`, id = "master")
)

outfile = normpath(@__DIR__,"results_$(target).json")
export_markdown(outfile,results)
