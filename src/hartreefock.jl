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

function rhf(basis, molecule::Molecule, maxiter = 20, convergence = 1e-6, maxdiis = 5)
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

    error_vectors = []
    Fock_list = []

    println("Starting SCF iterations...")

    for iteration = 0:maxiter
        Eold = Eel

        P = zeros(K, K)
        
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
        
        error = F * D * S - S * D * F
        push!(error_vectors, error)
        push!(Fock_list, F)

        if length(error_vectors) > maxdiis
            popfirst!(error_vectors)
            popfirst!(Fock_list)
        end

        if length(error_vectors) > 1
            dim_B = length(error_vectors) + 1
            B = zeros(dim_B, dim_B)
            B[dim_B, dim_B] = 0.0

            for i in 1:length(error_vectors), j in 1:length(error_vectors)
                B[i, j] = sum(error_vectors[i] .* error_vectors[j])
                B[dim_B, i] = -1.0
                B[i, dim_B] = -1.0
            end

            rhs = zeros(dim_B)
            rhs[dim_B] = -1.0
            coeffs = B \ rhs

            F = sum(coeffs[i] * Fock_list[i] for i in 1:length(error_vectors))
        end

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