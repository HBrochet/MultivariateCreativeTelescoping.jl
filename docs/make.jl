using Documenter, MultivariateCreativeTelescoping

DocMeta.setdocmeta!(MultivariateCreativeTelescoping, :DocTestSetup, :(using MultivariateCreativeTelescoping); recursive=true)

makedocs(sitename="MultivariateCreativeTelescoping.jl",
        pages = [
            "Presentation" => "index.md",
            "Quick start" => "QuickStart.md"
        ])

deploydocs(
    repo = "github.com/HBrochet/MultivariateCreativeTelescoping.jl.git",
)