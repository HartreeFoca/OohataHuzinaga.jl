struct Operator<:Function
    ô
end

function Kxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ)
    K = αⱼ * (2 * (ℓⱼ + mⱼ + nⱼ) + 3) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ)
    K -= (2 * (αⱼ^2)) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ + 2, mᵢ, mⱼ, nᵢ, nⱼ)
    K -= (2 * (αⱼ^2)) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ + 2, nᵢ, nⱼ)
    K -= (2 * (αⱼ^2)) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ + 2)
    K -= (1 / 2) * (ℓⱼ * (ℓⱼ - 1)) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ - 2, mᵢ, mⱼ, nᵢ, nⱼ)
    K -= (1 / 2) * (mⱼ * (mⱼ - 1)) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ - 2, nᵢ, nⱼ)
    K -= (1 / 2) * (nⱼ * (nⱼ - 1)) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ - 2)
    return K
end

function oei(basis, molecule::Molecule, integral::Operator)
    n = length(basis)
    O = zeros(n, n)
    for c in CartesianIndices(O)
        i, j = c[1], c[2]
        basisᵢ, basisⱼ = basis[i], basis[j]
        Rᵢ, Rⱼ = basisᵢ.R, basisⱼ.R
        dist = distance(Rᵢ, Rⱼ)
        m, p = basisᵢ.size, basisⱼ.size
        ℓᵢ, mᵢ, nᵢ = basisᵢ.ℓ, basisᵢ.m, basisᵢ.n
        ℓⱼ, mⱼ, nⱼ = basisⱼ.ℓ, basisⱼ.m, basisⱼ.n
        for e in CartesianIndices((m, p))
            k, l = e[1], e[2]
            αᵢ = basisᵢ.α[k]
            αⱼ = basisⱼ.α[l]
            dᵢ, dⱼ = basisᵢ.d[k], basisⱼ.d[l]
            Nᵢ, Nⱼ = basisᵢ.N[k], basisⱼ.N[l]
            O[i, j] += (
                exp(-αᵢ * αⱼ * dist / (αᵢ + αⱼ)) *
                Nᵢ * Nⱼ * dᵢ * dⱼ *
                integral.ô(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ)
            )
        end
    end
    return O
end

function kinetic(basis, molecule::Molecule)
    n = length(basis)
    T = zeros(n, n)

    for c in CartesianIndices(T)
        i, j = c[1], c[2]

        basisᵢ, basisⱼ = basis[i], basis[j]
        Rᵢ, Rⱼ = basisᵢ.R, basisⱼ.R

        dist = distance(Rᵢ, Rⱼ)

        m, p = basisᵢ.size, basisⱼ.size

        ℓᵢ, mᵢ, nᵢ = basisᵢ.ℓ, basisᵢ.m, basisᵢ.n
        ℓⱼ, mⱼ, nⱼ = basisⱼ.ℓ, basisⱼ.m, basisⱼ.n

        for e in CartesianIndices((m, p))
            k, l = e[1], e[2]

            αᵢ = basisᵢ.α[k]
            αⱼ = basisⱼ.α[l]

            dᵢ, dⱼ = basisᵢ.d[k], basisⱼ.d[l]
            Nᵢ, Nⱼ = basisᵢ.N[k], basisⱼ.N[l]

            T[i, j] += (
                exp(-αᵢ * αⱼ * dist / (αᵢ + αⱼ)) *
                Nᵢ * Nⱼ * dᵢ * dⱼ *
                Kxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ)
            )
        end
    end

    return T
end