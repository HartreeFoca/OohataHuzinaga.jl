@time @testset "kinetic.jl" begin
    methane = molecule("data/methane.xyz")
    basis = parsebasis(methane, "sto-3g")
    @test isapprox(
        kinetic(basis, methane),
        [
             15.8911      -0.08589      5.23466e-22   2.73557e-21  -4.10859e-21   0.093943  0.0939292   0.0939287   0.093935;
            -0.08589       0.47225      0.0           0.0           5.97539e-22   0.337973  0.337961    0.33796     0.337966;
            -4.96639e-21   0.0          1.47773       0.0           0.0          -0.104397  0.401344    0.186147   -0.483147;
            -3.00562e-21   0.0          0.0           1.47773       0.0           0.459873  0.130077   -0.43217    -0.157803;
             5.78625e-21  -4.4878e-22   0.0           0.0           1.47773      -0.328124  0.389921   -0.329566    0.267782;
             0.093943      0.337973    -0.104397      0.459873     -0.328124      0.760032  0.131392    0.131385    0.131408;
             0.0939292     0.337961     0.401344      0.130077      0.389921      0.131392  0.760032    0.131388    0.131387;
             0.0939287     0.33796      0.186147     -0.43217      -0.329566      0.131385  0.131388    0.760032    0.131401;
             0.093935      0.337966    -0.483147     -0.157803      0.267782      0.131408  0.131387    0.131401    0.760032
        ],
        atol = 1e-4,
    )
end
