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

function computeenergy(basis, molecule::Molecule, maxiter = 20, convergence = 1e-6)
    S = overlap(basis)
    T = kinetic(basis, molecule)
    V = attraction(basis, molecule)
    G = repulsion(basis, molecule)

    K = length(basis)

    Hcore = T .+ V

    D = zeros(K, K)
    P = zeros(K, K)

    X = sqrt(inv(S))
    println(X)

    Eel = 0.0

    for iteration = 0:maxiter
        Eold = Eel
        for n = 1:K
            for m = 1:K
                P[m, n] = 0.0
                for ℓ = 1:K
                    for s = 1:K
                        P[m, n] += D[ℓ, s] * (G[m, n, s, ℓ] - 0.5 * G[m, ℓ, s, n])
                    end
                end
            end
        end

        F = Hcore + P         
        Fp = X * F * X  

        eigen_decomp = eigen(Fp)   
        e = eigen_decomp.values

        Cp = eigen_decomp.vectors
        C = X * Cp

        N = electroncount(molecule)

        for n = 1:K
            for m = 1:K
                D[m, n] = 0.0
                for a = 1:trunc(Int64, N / 2)
                    D[m, n] += 2 * (C[m, a] * C[n, a])
                end
            end
        end

        Eel = 0.0

        for m = 1:K
            for n = 1:K
                Eel += 0.5 * D[n, m] * (Hcore[m, n] + F[m, n])
            end
        end

        if (abs(Eel - Eold) < convergence) && (iteration > 0)
            break
        end

        Vnn = nuclearrepulsion(molecule)

        println(Eel + Vnn)
    end
end
