
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
  "Gridap.ODEs" => "ODEs.md",
  "Gridap.Adaptivity" => "Adaptivity.md",
  "Developper notes" => Any[
    "dev-notes/block-assemblers.md",
  ],
]

makedocs(
  sitename = "Gridap.jl",
  format = Documenter.HTML(
    size_threshold=nothing
  ),
  modules = [Gridap],
  pages = pages,
  doctest = false,
  warnonly = [:cross_references,:missing_docs],
  checkdocs = :exports,
)

deploydocs(
  repo = "github.com/gridap/Gridap.jl.git",
)
