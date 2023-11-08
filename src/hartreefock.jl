struct Results
    energy::Float64
    orbitals
    density
end

function electroncount(molecule::Molecule)
    N = 0

    for atom in getatoms(molecule)
        N += atom.number
    end

    return N
end

function nuclearrepulsion(molecule::Molecule)
    Vnn = 0.0
    atoms = getatoms(molecule)

    for i in 1:length(atoms), j in 1:length(atoms)
        if j > i
            num = atoms[i].number * atoms[j].number
            den = sqrt(distance(atoms[i].coords, atoms[j].coords))
            Vnn += num / den
        end
    end

    return Vnn
end

function rhf(basis, molecule::Molecule, maxiter = 20, convergence = 1e-6)
    S = overlap(basis)
    println("Overlap is done!")

    T = kinetic(basis, molecule)
    println("Kinetic is done!")

    V = attraction(basis, molecule)
    println("Attraction is done!")

    G = repulsion(basis, molecule)
    println("Repulsion is done!")

    K = length(basis)

    Hcore = T .+ V
    println("HCore is done!")

    D = zeros(K, K)

    X = sqrt(inv(S))

    Eel = 0.0

    println("Starting SCF iterations...")

    for iteration = 0:maxiter
        Eold = Eel

        P = zeros(K, K)

        for m in 1:K, n in 1:K
            P[m, n] = sum(D[ℓ, s] * (G[m, n, s, ℓ] - 0.5 * G[m, ℓ, s, n]) for ℓ in 1:K, s in 1:K)
        end

        F = Hcore + P
        
        Fp = X * F * X  

        eigen_decomp = eigen(Fp)   
        e = eigen_decomp.values

        Cp = eigen_decomp.vectors
        C = X * Cp

        N = electroncount(molecule)

        D = [2 * sum(C[m, a] * C[n, a] for a in 1:trunc(Int64, N / 2)) for m in 1:K, n in 1:K]

        Eel = 0.5 * sum(D[n, m] * (Hcore[m, n] + F[m, n]) for m in 1:K, n in 1:K)

        if (abs(Eel - Eold) < convergence) && (iteration > 0)
            break
        end

        Vnn = nuclearrepulsion(molecule)

        println(Eel + Vnn)
        println(Vnn)
    end

    return Results(Eel, D)
end