using Documenter, TermoQuimica

DocMeta.setdocmeta!(TermoQuimica, :DocTestSetup, :(using TermoQuimica); recursive=true)


makedocs(
    repo="https://github.com/EmilioAlvizo/TermoQuimica.jl/blob/{commit}{path}#{line}",
    sitename="TermoQuimica.jl",
    modules = [TermoQuimica],
    authors  = "Emilio Alvizo Velázquez",

    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical="https://EmilioAlvizo.github.io/TermoQuimica.jl",
        assets = String[]),

    pages = [
        "Inicio" => "index.md",
        "Background" => "background.md",
        "Guia de usuario" => "userguide.md",
        "API" => "api.md",
        "Servidor" => "servidor.md",
    ])

deploydocs(
    repo="github.com/EmilioAlvizo/TermoQuimica.jl",
)