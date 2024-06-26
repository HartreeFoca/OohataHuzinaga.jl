@time @testset "overlap.jl" begin
    methane = BasisSets.molecule("data/methane.xyz")
    basis = BasisSets.parsebasis(methane, "sto-3g")
    @test isapprox(
        overlap(basis),
        [
            1.0          0.248362      7.39808e-23   9.5106e-22  -4.61187e-22   0.193855   0.193845   0.193845   0.193849;
            0.248362     1.0           0.0           0.0          1.01855e-21   0.794111   0.794097   0.794096   0.794103;
            7.39808e-23  0.0           1.0           0.0          0.0          -0.0939537  0.361208   0.167531  -0.434825;
            9.5106e-22   0.0           0.0           1.0          0.0           0.41387    0.117069  -0.388952  -0.14202;
           -4.61187e-22  1.01855e-21   0.0           0.0          1.0          -0.295301   0.350927  -0.296608   0.240999;
            0.193855     0.794111     -0.0939537     0.41387     -0.295301      1.0        0.5271     0.527089   0.527123;
            0.193845     0.794097      0.361208      0.117069     0.350927      0.5271     1.0        0.527093   0.527093;
            0.193845     0.794096      0.167531     -0.388952    -0.296608      0.527089   0.527093   1.0        0.527113;
            0.193849     0.794103     -0.434825     -0.14202      0.240999      0.527123   0.527093   0.527113   1.0
        ],
        atol = 1e-4,
    )
end
