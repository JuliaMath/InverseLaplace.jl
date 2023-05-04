using Documenter, InverseLaplace

makedocs(;
         modules=[InverseLaplace],
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
#         pages=["Home" => "index.md"],
         pages=["index.md"],
         repo="https://github.com/JuliaMath/InverseLaplace.jl/blob/{commit}{path}#L{line}",
         sitename="InverseLaplace.jl",
         authors="John Lapeyre",
)

deploydocs(; repo="github.com/JuliaMath/InverseLaplace.jl", push_preview=true)

# using Documenter, InverseLaplace
# makedocs(
#     format = :html,
#     sitename = "InverseLaplace.jl",
#     modules = [InverseLaplace],
#     pages = [
#         "index.md"
#     ]
# )

# deploydocs(
#     repo = "github.com/JuliaMath/InverseLaplace.jl.git",
#     target = "build",
#     julia  = "1.8",
#     osname = "linux",
#     deps = nothing,
#     make = nothing
# )
