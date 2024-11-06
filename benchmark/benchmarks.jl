using BenchmarkTools
using PkgBenchmark
using Gridap

macro include_bm(SUITE,name)
  quote
    include("bm/$($name).jl")
    SUITE["$($name)"] = $(Symbol(name)).SUITE
  end
end

const SUITE = BenchmarkGroup()

@include_bm SUITE "bm_assembly"
