using OohataHuzinaga
using Test

@testset "OohataHuzinaga.jl" begin
    include("molecule.jl")
    include("basis.jl")
    include("doublefactorial.jl")
    include("normalization.jl")
    include("overlap.jl")
    include("kinetic.jl")
    include("attraction.jl")
    include("repulsion.jl")
end
