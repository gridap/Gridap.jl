using Documenter
using Gridap

pages = [
  "Home" => "index.md",
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
    sitename = "Gridap",
    format = Documenter.HTML(),
    modules = [Gridap],
    pages = pages
)

deploydocs(
    repo = "github.com/gridap/Gridap.jl.git",
)

