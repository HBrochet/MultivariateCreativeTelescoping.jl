using Documenter, MultivariateCreativeTelescoping

DocMeta.setdocmeta!(MultivariateCreativeTelescoping, :DocTestSetup, :(using MultivariateCreativeTelescoping); recursive=true)

makedocs(sitename="MultivariateCreativeTelescoping.jl",
        modules=[MultivariateCreativeTelescoping],
        format = Documenter.HTML(),
        pages = [
            "Package presentation" => "index.md",
            "Quick start" => "QuickStart.md",
            "Ore algebra" => "OreAlgebra.md",
            "Multivariate Creative Telescoping" =>"MCT.md",
            "Other functionalities" => "OtherFunctionalities.md",

        ],
        checkdocs=:none
        )

deploydocs(
    repo = "github.com/HBrochet/MultivariateCreativeTelescoping.jl",
    target = "build",
    branch = "gh-pages",
)