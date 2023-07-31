using CrystallographyBase
using Documenter

DocMeta.setdocmeta!(CrystallographyBase, :DocTestSetup, :(using CrystallographyBase); recursive=true)

makedocs(;
    modules=[CrystallographyBase],
    authors="singularitti <singularitti@outlook.com> and contributors",
    repo="https://github.com/MineralsCloud/CrystallographyBase.jl/blob/{commit}{path}#{line}",
    sitename="CrystallographyBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/CrystallographyBase.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Installation Guide" => "man/installation.md",
            "Definitions and conventions" => "man/definitions.md",
            "Examples" => "man/examples.md",
            "Troubleshooting" => "man/troubleshooting.md",
        ],
        "Reference" => Any[
            "Public API" => "lib/public.md",
            # "Internals" => map(
            #     s -> "lib/internals/$(s)",
            #     sort(readdir(joinpath(@__DIR__, "src/lib/internals")))
            # ),
        ],
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style-guide.md",
            "Design Principles" => "developers/design-principles.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/CrystallographyBase.jl",
    devbranch="main",
)
