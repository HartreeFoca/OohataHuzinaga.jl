using BasisSets
using OohataHuzinaga
using Test
@testset "OohataHuzinaga.jl" begin
    include("overlap.jl")
    include("kinetic.jl")
    include("attraction.jl")
end
