using OohataHuzinaga
using Test
using BasisSets

@testset "OohataHuzinaga.jl" begin
    include("doublefactorial.jl")
    include("normalization.jl")
    include("overlap.jl")
    include("kinetic.jl")
    include("attraction.jl")
    include("repulsion.jl")
end
