push!(LOAD_PATH,"../src/")
using Documenter 
using Gridap
makedocs(
         sitename = "Gridap",
         format = Documenter.HTML(),
         modules = [Gridap]
        )

deploydocs(
           repo = "github.com/gridap/Gridap.jl.git"
          )

