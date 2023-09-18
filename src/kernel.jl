"""
This function calculates the integral kernel
"""
function kernel(basis, operator)
    n = length(basis) 
    S = zeros(n, n)

    for i in 1:n, j in 1:n
        basisᵢ = basis[i]
        basisⱼ = basis[j]
        
        Rᵢ = basisᵢ.R
        Rⱼ = basisⱼ.R

        dist = distance(Rᵢ, Rⱼ)

        m = length(basisᵢ.α)
        p = length(basisⱼ.α)

        @views for k in 1:m, l in 1:p
            αᵢ = basisᵢ.α[k]
            αⱼ = basisⱼ.α[l]

            dᵢ = basisᵢ.d[k]
            dⱼ = basisⱼ.d[l]

            ℓᵢ, mᵢ, nᵢ = basisᵢ.ℓ, basisᵢ.m, basisᵢ.n
            ℓⱼ, mⱼ, nⱼ = basisⱼ.ℓ, basisⱼ.m, basisⱼ.n

            operator(αᵢ, αⱼ, Rᵢ, Rⱼ, dist, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ)

            S[i, j] += (
                exp(-αᵢ * αⱼ * dist / (αᵢ + αⱼ)) *
                normalization(αᵢ, ℓᵢ, mᵢ, nᵢ) *
                normalization(αⱼ, ℓⱼ, mⱼ, nⱼ) *
                dᵢ *
                dⱼ *
                Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ)
            )
        end
    end

    return S
end