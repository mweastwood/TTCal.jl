using Documenter, TTCal
makedocs(
    format = :html,
    sitename = "TTCal.jl",
    authors = "Michael Eastwood",
    linkcheck = true,
    html_prettyurls = !("local" in ARGS),
    pages = [
        "Home" => "index.md",
        "Internal API" => "internal-api.md",
    ]
)

deploydocs(
    repo = "github.com/mweastwood/TTCal.jl.git",
    julia = "0.6",
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing
)

