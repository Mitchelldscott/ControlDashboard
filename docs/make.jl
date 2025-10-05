using Documenter
using ControlDashboard 

# 1. Define the build settings for the documentation.
# The `modules` keyword tells Documenter which modules to analyze for docstrings.
# The `pages` keyword defines the structure/TOC of your documentation site.
makedocs(
    modules = [ControlDashboard, ControlPanel],
    sitename = "ControlDashboard.jl Documentation",
    authors = "Mitchell D Scott",
    clean = true, 
    format = Documenter.HTML(
        # Set prettyurls=true for nice URLs when hosted (e.g., on GitHub Pages)
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://your_github_username.github.io/ControlDashboard.jl/stable/", # Update this URL
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Introduction" => "man/introduction.md",
            "Examples" => "man/examples.md",
        ],
        "API Reference" => "api.md",
        # Uncomment this line if you want to include all unlisted docstrings
        # "Index" => "function_index.md",
    ]
)

# 2. Deploy the documentation to GitHub Pages.
# This step is conditional on running within a Continuous Integration environment.
deploydocs(
    repo = "github.com/your_github_username/ControlDashboard.jl.git", # Update with your repo
    devbranch = "master", # The branch where development happens (should match your CI.yml 'push' branch)
    push_preview = true, # Allows for deploying docs for pull requests
)

# You can also add Pkg.activate(dirname(@__DIR__)) if you need to ensure the docs
# environment is active, though the CI script should handle this with '--project=.'
