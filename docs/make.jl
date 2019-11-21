
using Documenter
using Gridap

pages = [
  "Home" => "index.md",
  "Getting Started" => "getting-started.md",
  "Gridap" => "Gridap.md",
  "Gridap.Helpers" => "Helpers.md",
  "Gridap.Inference" => "Inference.md",
  "Gridap.TensorValues" => "TensorValues.md",
  "Gridap.Arrays" => "Arrays.md",
  "Gridap.Fields" => "Fields.md",
  "Gridap.Polynomials" => "Polynomials.md",
  "Gridap.ReferenceFEs" => "ReferenceFEs.md",
 ]

makedocs(
    sitename = "Gridap.jl",
    format = Documenter.HTML(),
    modules = [Gridap],
    pages = pages
)

deploydocs(
    repo = "github.com/gridap/Gridap.jl.git",
)

#using Documenter 
#using Gridap
#
#makedocs(
#  sitename = "Gridap.jl",
#  format = Documenter.HTML(),
#  modules = [Gridap],
#  pages = [
#    "Home" => "index.md",
#    "Getting Started" => "pages/getting-started.md",
#    "Manual" => "pages/manual.md"]
# )
#
#deploydocs(repo = "github.com/gridap/Gridap.jl.git")
#
