using Documenter
include("../src/FSI.jl")
using FSI

makedocs(
    format = :html,
    sitename = "FSI.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => ["manual/problemsetup.md",
                    # "manual/timemarching.md",
                    # "manual/saddlepointsystem.md",
                     ]
    ],
    assets = ["assets/custom.css"],
    strict = true
)


# if "DOCUMENTER_KEY" in keys(ENV)
deploydocs(
 repo = "github.com/ruizhi92/FSI.jl.git",
 target = "build",
 branch = "gh-pages",
 deps = nothing,
 make = nothing,
 julia = "1.3"
)
# end
