push!(LOAD_PATH,"../src/")
using Documenter 
using Numa
makedocs(
         sitename = "Numa",
         format = Documenter.HTML(),
         modules = [Numa]
        )

deploydocs(
           repo = "github.com/lssc-team/Numa.jl.git"
          )

