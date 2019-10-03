using Documenter 
using Gridap

makedocs(
  sitename = "Gridap.jl",
  format = Documenter.HTML(),
  modules = [Gridap],
  pages = [
    "Home" => "index.md",
    "Getting Started" => "pages/getting-started.md",
    "Manual" => "pages/manual.md"]
 )

deploydocs(repo = "github.com/gridap/Gridap.jl.git")

