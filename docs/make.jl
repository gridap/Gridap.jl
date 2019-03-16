push!(LOAD_PATH,"../src/")
using Documenter 
using Numa
makedocs(
         sitename = "Numa",
         format = Documenter.HTML(),
         modules = [Numa]
        )

deploydocs(
           repo = "github.com/LSSC-team/Numa.jl.git"
          )

