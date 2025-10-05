using Pkg; Pkg.develop(path = "."); Pkg.instantiate()

using Documenter
using ControlDashboard 

# 1. Define the build settings for the documentation.
# The `modules` keyword tells Documenter which modules to analyze for docstrings.
# The `pages` keyword defines the structure/TOC of your documentation site.
makedocs(
    modules = [ControlDashboard],
    sitename = "ControlDashboard.jl Documentation",
    authors = "Mitchell D Scott",
    clean = true, 
    format = Documenter.HTML(
        # Set prettyurls=true for nice URLs when hosted (e.g., on GitHub Pages)
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://MitchellDScott.github.io/ControlDashboard/stable/",
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API Reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/MitchellDScott/ControlDashboard.git",
    devbranch = "master",
    push_preview = true,
)