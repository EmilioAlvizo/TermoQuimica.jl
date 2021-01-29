push!(LOAD_PATH,"../src/")

using Documenter, TermoQuimica

makedocs(
    sitename="Documentación",
    format = Documenter.HTML(prettyurls = false,assets = ["assets/custom.css"]),
    modules = [TermoQuimica],
    authors  = "Emilio Alvizo Velázquez",
    pages = [
        "Inicio" => "index.md",
        "Background" => "background.md",
        "Guia de usuario" => "userguide.md",
        "API" => "api.md",
        "Servidor" => "servidor.md",
    ])