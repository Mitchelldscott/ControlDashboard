using Documenter, ControlDashboard

makedocs(
    sitename = "ControlDashboard.jl",
    modules = [ControlDashboard],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/MitchellDScott/ControlDashboard.jl.git",
)