using Documenter, InverseLaplace

makedocs(
    modules = InverseLaplace,
    clean   = false,
)

# deploydocs(
#     deps = Deps.pip("pygments", "mkdocs", "mkdocs-material"),
#     repo = "github.com/MichaelHatherly/PrivateModules.jl.git",
# )
