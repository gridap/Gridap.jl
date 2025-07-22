
using Documenter
using Gridap

pages = [
  "Home" => "index.md",
  "Getting Started" => "getting-started.md",
  "Gridap at a glance" => "Gridap.md",
  "Ecosystem" => "ecosystem.md",
  "Modules" => [
    "Helpers" => "Helpers.md",
    "Io" => "Io.md",
    "Algebra" => "Algebra.md",
    "Arrays" => "Arrays.md",
    "TensorValues" => "TensorValues.md",
    "Fields" => "Fields.md",
    "Polynomials" => "Polynomials.md",
    "ReferenceFEs" => "ReferenceFEs.md",
    "Geometry" => "Geometry.md",
    "CellData" => "CellData.md",
    "Visualization" => "Visualization.md",
    "FESpaces" => "FESpaces.md",
    "MultiField" => "MultiField.md",
    "ODEs" => "ODEs.md",
    "Adaptivity" => "Adaptivity.md",
  ],
  "Developper notes" => Any[
    "dev-notes/block-assemblers.md",
    "dev-notes/pullbacks.md",
    "dev-notes/bernstein.md",
    "dev-notes/autodiff.md",
  ],
]

makedocs(
  sitename = "Gridap.jl",
  format = Documenter.HTML(
    size_threshold=nothing,
    size_threshold_warn=300 * 2^10 # 300KiB
  ),
  modules = [Gridap],
  pages = pages,
  doctest = false,
  warnonly = [:missing_docs], # ,:cross_references
  checkdocs = :exports,
)

deploydocs(
  repo = "github.com/gridap/Gridap.jl.git",
)
