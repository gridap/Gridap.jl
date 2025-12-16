
using Documenter
using Gridap
using TikzPictures

pages = [
  "Home" => "index.md",
  "Getting Started" => "getting-started.md",
  "Gridap at a glance" => "overview.md",
  "Ecosystem" => "ecosystem.md",
  "Modules" => [
    "Helpers" => "modules/Helpers.md",
    "Io" => "modules/Io.md",
    "Algebra" => "modules/Algebra.md",
    "Arrays" => "modules/Arrays.md",
    "TensorValues" => "modules/TensorValues.md",
    "Fields" => "modules/Fields.md",
    "Polynomials" => "modules/Polynomials.md",
    "ReferenceFEs" => "modules/ReferenceFEs.md",
    "Geometry" => "modules/Geometry.md",
    "CellData" => "modules/CellData.md",
    "Visualization" => "modules/Visualization.md",
    "FESpaces" => "modules/FESpaces.md",
    "MultiField" => "modules/MultiField.md",
    "ODEs" => "modules/ODEs.md",
    "Adaptivity" => "modules/Adaptivity.md",
  ],
  "Extensions" => [
    "TikzPictures" => "extensions/TikzPictures.md",
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
    size_threshold_warn=500 * 2^10 # 300KiB
  ),
  modules = [Gridap],
  pages = pages,
  doctest = false,
  warnonly = [:missing_docs,:cross_references], # ,:cross_references
  checkdocs = :exports,
)

deploydocs(
  repo = "github.com/gridap/Gridap.jl.git",
)
