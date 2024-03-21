using Documenter, MultivariateCreativeTelescoping

DocMeta.setdocmeta!(MultivariateCreativeTelescoping, :DocTestSetup, :(using MultivariateCreativeTelescoping); recursive=true)

makedocs(sitename="MultivariateCreativeTelescoping.jl",
        modules=[MultivariateCreativeTelescoping],
        pages = [
            "Presentation" => "index.md",
            "Quick start" => "QuickStart.md"
        ],
        checkdocs=:none
        )

deploydocs(
    repo = "github.com/HBrochet/MultivariateCreativeTelescoping.jl",
    target = "build",
    branch = "gh-pages",
)