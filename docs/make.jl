push!(LOAD_PATH, "../src/")

using Ohata
using Documenter

makedocs(
    sitename = "Ohata.jl",
    modules = [Ohata],
    pages = [
        "Home" => "index.md",
        "API" => ["Input" => "input.md", "Basis Sets" => "basis.md"],
    ],
)

deploydocs(; repo = "github.com/HartreeFoca/Ohata.jl.git")
