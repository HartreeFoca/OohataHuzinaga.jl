@time @testset "repulsion.jl" begin
    hydrogen = BesisSets.molecule("data/h2.xyz")
    basis = BasisSets.parsebasis(hydrogen, "sto-3g")
    @test isapprox(
        repulsion(basis, hydrogen),
        [
            0.7746059298062673 0.6482052325853882
            0.6482052325853882 0.6994796252934075;;;
            0.6482052325853882 0.5707519808185033
            0.5707519808185033 0.6482052325853881;;;;
            0.6482052325853882 0.5707519808185033
            0.5707519808185033 0.6482052325853881;;;
            0.6994796252934075 0.6482052325853882
            0.6482052325853883 0.7746059298062673
        ],
        atol = 1e-4,
    )
end
