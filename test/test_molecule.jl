@time @testset "molecule.jl" begin
    methane = molecule("methane.xyz")
    @test methane.atoms  == ["C", "H", "H", "H", "H"]
    @test methane.coords == [0.00001021434087  0.00001532972083 -0.00001493500137;
                            -0.19951695340554  0.87894179053067 -0.62713882127936;
                             0.76712229809243  0.24863902907755  0.74526241504934;
                             0.35580334399536 -0.82601803138729 -0.62993342769733;
                            -0.92343260142312 -0.30159515034176  0.51179839372872]
end