@time @testset begin
    methane = BasisSets.molecule("data/methane.xyz")
    sto3gBasis = BasisSets.parsebasis(methane, "sto-3g")
    @test sto3gBasis[1].R == [0.00001021434087, 0.00001532972083, -0.00001493500137]
    @test sto3gBasis[1].α == [3.53051220, 13.0450960, 71.6168370]
    @test sto3gBasis[1].d == [0.44463454, 0.53532814, 0.15432897]
end
