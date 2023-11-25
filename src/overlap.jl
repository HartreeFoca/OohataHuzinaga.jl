function Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ)
    Rₚ = gaussianproduct(αᵢ, Rᵢ, αⱼ, Rⱼ, αᵢ + αⱼ)

    Sx = sᵢ(ℓᵢ, ℓⱼ, αᵢ + αⱼ, Rᵢ[1], Rⱼ[1], Rₚ[1])
    Sy = sᵢ(mᵢ, mⱼ, αᵢ + αⱼ, Rᵢ[2], Rⱼ[2], Rₚ[2])
    Sz = sᵢ(nᵢ, nⱼ, αᵢ + αⱼ, Rᵢ[3], Rⱼ[3], Rₚ[3])

    return Sx * Sy * Sz
end

function overlap(basis)
    n = length(basis)
    S = zeros(n, n)

    for c in CartesianIndices(S)
        i, j = c[1], c[2]

        basisᵢ, basisⱼ = basis[i], basis[j]
        Rᵢ, Rⱼ = basisᵢ.R, basisⱼ.R
        m, p = basisᵢ.size, basisⱼ.size

        ℓᵢ, mᵢ, nᵢ = basisᵢ.ℓ, basisᵢ.m, basisᵢ.n
        ℓⱼ, mⱼ, nⱼ = basisⱼ.ℓ, basisⱼ.m, basisⱼ.n

        dist = distance(Rᵢ, Rⱼ)
        
        for e in CartesianIndices((m, p))
            k, l = e[1], e[2]

            αᵢ, αⱼ = basisᵢ.α[k], basisⱼ.α[l]
            dᵢ, dⱼ = basisᵢ.d[k], basisⱼ.d[l]

            Nᵢ, Nⱼ = basisᵢ.N[k], basisⱼ.N[l]

            S[i, j] += (
                exp(-αᵢ * αⱼ * dist / (αᵢ + αⱼ)) *
                Nᵢ * Nⱼ * dᵢ * dⱼ * 
                Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ)
            )
        end
    end

    return S
end
