using Documenter, MendelGeneDropping

ENV["DOCUMENTER_DEBUG"] = "true"

makedocs(
    format = :html,
    sitename = "MendelGeneDropping",
    modules = [MendelGeneDropping],
    authors = "Jeanette Papp",
    clean = true,
    debug = true,
    pages = [
        "index.md"
    ]
)

deploydocs(
    repo   = "github.com/OpenMendel/MendelGeneDropping.jl.git",
    target = "build",
    julia  = "1.0",
    deps   = nothing,
    make   = nothing
)
