
using Documenter
using Gridap

pages = [
  "Home" => "index.md",
  "Getting Started" => "getting-started.md",
  "Gridap" => "Gridap.md",
  "Gridap.Helpers" => "Helpers.md",
  "Gridap.Io" => "Io.md",
  "Gridap.Algebra" => "Algebra.md",
  "Gridap.Arrays" => "Arrays.md",
  "Gridap.TensorValues" => "TensorValues.md",
  "Gridap.Fields" => "Fields.md",
  "Gridap.Polynomials" => "Polynomials.md",
  "Gridap.ReferenceFEs" => "ReferenceFEs.md",
  "Gridap.Geometry" => "Geometry.md",
  "Gridap.CellData" => "CellData.md",
  "Gridap.Visualization" => "Visualization.md",
  "Gridap.FESpaces" => "FESpaces.md",
  "Gridap.MultiField" => "MultiField.md",
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

