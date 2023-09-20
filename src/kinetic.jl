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

function kinetic(basis, molecule::Molecule)
    n = length(basis)
    T = zeros(n, n)

    for i in 1:n, j in 1:n
        basisᵢ = basis[i]
        basisⱼ = basis[j]
        
        Rᵢ = basisᵢ.R
        Rⱼ = basisⱼ.R

        dist = distance(Rᵢ, Rⱼ)

        m = length(basisᵢ.α)
        p = length(basisⱼ.α)

        for k in 1:m, l in 1:p
            αᵢ = basisᵢ.α[k]
            αⱼ = basisⱼ.α[l]

            dᵢ = basisᵢ.d[k]
            dⱼ = basisⱼ.d[l]

            ℓᵢ, mᵢ, nᵢ = basisᵢ.ℓ, basisᵢ.m, basisᵢ.n
            ℓⱼ, mⱼ, nⱼ = basisⱼ.ℓ, basisⱼ.m, basisⱼ.n

            T[i, j] += (
                exp(-αᵢ * αⱼ * dist / (αᵢ + αⱼ)) *
                normalization(αᵢ, ℓᵢ, mᵢ, nᵢ) *
                normalization(αⱼ, ℓⱼ, mⱼ, nⱼ) *
                dᵢ *
                dⱼ *
                Kxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ)
            )
        end
    end

    return T
end