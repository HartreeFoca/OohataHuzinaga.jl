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
                                        Vz = vᵢ(
                                            n,
                                            t,
                                            k,
                                            nᵢ,
                                            nⱼ,
                                            Rᵢ[3],
                                            Rⱼ[3],
                                            Rₖ[3],
                                            Rₚ[3],
                                            γ,
                                        )

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

    for i in 1:n, j in 1:n
        basisᵢ = basis[i]
        basisⱼ = basis[j]
        
        Rᵢ = basisᵢ.R
        Rⱼ = basisⱼ.R

        m = length(basisᵢ.α)
        p = length(basisⱼ.α)

        for natom in 1:natoms
            Rₖ = molecule.coords[natom, :]

            el = 0
            for k in 1:m, l in 1:p
                αᵢ = basisᵢ.α[k]
                αⱼ = basisⱼ.α[l]
    
                dᵢ = basisᵢ.d[k]
                dⱼ = basisⱼ.d[l]
    
                ℓᵢ, mᵢ, nᵢ = basisᵢ.ℓ, basisᵢ.m, basisᵢ.n
                ℓⱼ, mⱼ, nⱼ = basisⱼ.ℓ, basisⱼ.m, basisⱼ.n

                Nᵢ = basisᵢ.N[k]
                Nⱼ = basisⱼ.N[l]

                el += Nᵢ * Nⱼ * dᵢ * dⱼ * Vxyz(ℓᵢ, mᵢ, nᵢ, ℓⱼ, mⱼ, nⱼ, αᵢ, αⱼ, Rᵢ, Rⱼ, Rₖ, atomicnumbers[natom])
            end
            V[i, j, natom] = el
        end
    end
    
    return sum(V, dims=3)[:,:,1]
end

function reducedimensions(matrix)
    r = zeros(eltype(matrix), axes(matrix, 1), axes(matrix, 2))
    sum!(r, matrix)
end