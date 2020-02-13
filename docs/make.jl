
using Documenter
using Gridap

pages = [
  "Home" => "index.md",
  "Getting Started" => "getting-started.md",
  "Gridap" => "Gridap.md",
  "Gridap.Helpers" => "Helpers.md",
  "Gridap.Inference" => "Inference.md",
  "Gridap.Algebra" => "Algebra.md",
  "Gridap.Io" => "Io.md",
  "Gridap.TensorValues" => "TensorValues.md",
  "Gridap.Arrays" => "Arrays.md",
  "Gridap.Fields" => "Fields.md",
  "Gridap.Polynomials" => "Polynomials.md",
  "Gridap.Integration" => "Integration.md",
  "Gridap.ReferenceFEs" => "ReferenceFEs.md",
  "Gridap.Geometry" => "Geometry.md",
  "Gridap.FESpaces" => "FESpaces.md",
  "Gridap.MultiField" => "MultiField.md",
  "Gridap.Visualization" => "Visualization.md",
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

