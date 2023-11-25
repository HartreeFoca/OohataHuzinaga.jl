function vᵢ(l, r, i, ℓᵢ, ℓⱼ, Aᵢ, Bᵢ, Cᵢ, Pᵢ, γ)
    ϵ = 1 / (4 * γ)

    vᵢ = (-1)^l
    vᵢ *= cₖ(l, ℓᵢ, ℓⱼ, Pᵢ - Aᵢ, Pᵢ - Bᵢ)
    vᵢ *= (-1)^i * factorial(l)
    vᵢ *= (Pᵢ - Cᵢ)^(l - (2 * r) - (2 * i)) * ϵ^(r + i)
    vᵢ /= factorial(r)
    vᵢ /= factorial(i)
    vᵢ /= factorial(l - (2 * r) - (2 * i))

    return vᵢ
end

function Vxyz(ℓᵢ, mᵢ, nᵢ, ℓⱼ, mⱼ, nⱼ, αᵢ, αⱼ, Rᵢ, Rⱼ, Rₖ, Z)
    γ = αᵢ + αⱼ

    Rₚ = gaussianproduct(αᵢ, Rᵢ, αⱼ, Rⱼ, γ)

    IJ = distance(Rᵢ, Rⱼ)
    PK = distance(Rₚ, Rₖ)

    Vxyz = 0
    for l = 0:(ℓᵢ+ℓⱼ)
        for r = 0:trunc(Int64, (l / 2))
            for i = 0:trunc(Int64, ((l - (2 * r)) / 2))
                Vx = vᵢ(l, r, i, ℓᵢ, ℓⱼ, Rᵢ[1], Rⱼ[1], Rₖ[1], Rₚ[1], γ)

                for m = 0:(mᵢ+mⱼ)
                    for s = 0:trunc(Int64, (m / 2))
                        for j = 0:trunc(Int64, ((m - (2 * s)) / 2))
                            Vy = vᵢ(m, s, j, mᵢ, mⱼ, Rᵢ[2], Rⱼ[2], Rₖ[2], Rₚ[2], γ)

                            for n = 0:(nᵢ+nⱼ)
                                for t = 0:trunc(Int64, (n / 2))
                                    for k = 0:trunc(Int64, ((n - (2 * t)) / 2))
                                        Vz = vᵢ(n, t, k, nᵢ, nⱼ, Rᵢ[3], Rⱼ[3], Rₖ[3], Rₚ[3], γ)

                                        ν = l + m + n - 2 * (r + s + t) - (i + j + k)
                                        F = boys(ν, (γ * abs(PK)))

                                        Vxyz += Vx * Vy * Vz * F
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    Vxyz *= (2 * π) / γ
    Vxyz *= exp(-αᵢ * αⱼ * abs(IJ) / γ)
    Vxyz *= -Z

    return Vxyz
end

function attraction(basis, molecule::Molecule)
    natoms = length(molecule.atoms)
    atomicnumbers = molecule.numbers

    n = length(basis)
    V = zeros(n, n, natoms)

    for c in CartesianIndices(V)
        i, j, natom = c[1], c[2], c[3]

        basisᵢ, basisⱼ = basis[i], basis[j]
        Rᵢ, Rⱼ = basisᵢ.R, basisⱼ.R

        ℓᵢ, mᵢ, nᵢ = basisᵢ.ℓ, basisᵢ.m, basisᵢ.n
        ℓⱼ, mⱼ, nⱼ = basisⱼ.ℓ, basisⱼ.m, basisⱼ.n

        Rₖ = molecule.coords[natom, :]
        m, p = basisᵢ.size, basisⱼ.size

        for e in CartesianIndices((m, p))
            k, l = e[1], e[2]

            αᵢ, αⱼ = basisᵢ.α[k], basisⱼ.α[l]
            dᵢ, dⱼ = basisᵢ.d[k], basisⱼ.d[l]
            Nᵢ, Nⱼ = basisᵢ.N[k], basisⱼ.N[l]

            V[i, j, natom] += Nᵢ * Nⱼ * dᵢ * dⱼ * Vxyz(ℓᵢ, mᵢ, nᵢ, ℓⱼ, mⱼ, nⱼ, αᵢ, αⱼ, Rᵢ, Rⱼ, Rₖ, atomicnumbers[natom])
        end
    end
    
    return sum(V, dims=3)[:,:,1]
end

function reducedimensions(matrix)
    r = zeros(eltype(matrix), axes(matrix, 1), axes(matrix, 2))
    sum!(r, matrix)
end