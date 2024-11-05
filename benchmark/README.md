# Benchmarking Suite

The following benchmarking suite uses `PkgBenchmark.jl` and `BenchmarkTools.jl` to compare the performance of different branches of Gridap. 

## Running the benchmarks

### Running as CI job

The benchmarks are setup as a manual Github Actions workflow in `.github/workflows/benchmark.yml`. To run the workflow, you will need administrator access to the Gridap repository, then follow instructions [here](https://docs.github.com/en/actions/managing-workflow-runs-and-deployments/managing-workflow-runs/manually-running-a-workflow).

The workflow will has two inputs: `target` and `base`, which are the branches/tags/commits you want to compare. The workflow will run the benchmarks on the `target` branch and compare them with the `base` branch (`master` by default).

You can also run the workflow using Github CLI and the following command:

```bash
gh workflow run benchmark.yml -f target=your_target_branch -f base=your_base_branch
```

### Running Locally

To run the benchmarks locally, you can have a look at the [documentation for `PkgBenchmark.jl`](https://juliaci.github.io/PkgBenchmark.jl/stable/run_benchmarks/).

Alternatively, you can run the CI script locally from a local copy of the repository. From the Gridap root directory, run the following commands:

```bash
# Instantiate Gridap and Gridap/benchmark
julia -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=benchmark/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
# Run the benchmarks
export BM_TARGET = "your target branch"
export BM_BASE = "your reference branch"
julia --project=benchmark/ --color=yes benchmark/run_benchmarks.jl
```

where `BM_TARGET` and `BM_BASE` are the branches/tags/commits you want to compare.

## Adding a new benchmark

To add a new benchmark suite `xyx`, create a new file `bm/bm_xyx.jl` that with the following structure:

```julia
module bm_xyx

using PkgBenchmark, BenchmarkTools

const SUITE = BenchmarkGroup()

[... Add your benchmarks here ...]

end # module
```

Then, add the following line to the `benchmarks.jl` file:

```julia
@include_bm "bm_xyz"
```

This should automatically include the new benchmarks into the global benchmarking suite.
