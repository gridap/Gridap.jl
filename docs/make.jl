using Documenter 
using Gridap

makedocs(
  sitename = "Gridap.jl",
  format = Documenter.HTML(),
  modules = [Gridap],
  pages = [
    "Home" => "index.md",
    "Getting Started" => "pages/getting-started.md",
    "Manual" => "pages/manual.md",
    "API"=>"pages/api.md",
    "Developer's Guide" => "pages/dev.md"]
 )

deploydocs(repo = "github.com/gridap/Gridap.jl.git")

