push!(LOAD_PATH, "../src/")

using OohataHuzinaga
using Documenter

makedocs(
    sitename = "OohataHuzinaga.jl",
    modules = [OohataHuzinaga],
    pages = [
        "Home" => "index.md",
        "API" => ["Input" => "input.md", "Basis Sets" => "basis.md"],
    ],
)

deploydocs(; repo = "github.com/HartreeFoca/OohataHuzinaga.jl.git")
