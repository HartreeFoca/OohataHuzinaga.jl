function doublefactorial(number)
    fact = foldl(Base.:*, range(number, 1, step=-2))

    return fact
end

function distance(Rᵢ, Rⱼ)
    d = (Rᵢ[1] - Rⱼ[1])^2 + (Rᵢ[2] - Rⱼ[2])^2 + (Rᵢ[3] - Rⱼ[3])^2

    return d
end

function normalization(α, ℓ, m, n)
    N = (4 * α)^(ℓ + m + n)
    N /=
        doublefactorial(2 * ℓ - 1) * doublefactorial(2 * m - 1) * doublefactorial(2 * n - 1)
    N *= ((2 * α) / π)^(3 / 2)
    N = sqrt(N)

    return N
end

function cₖ(j, l, m, A, B)
    coefficient = 0
    
    for k in 0:l, i in 0:m
        if (i + k == j)
            coefficient += binomial(l, k) * binomial(m, i) * A^(l - k) * B^(m - i)
        end
    end

    return coefficient
end

function sᵢ(ℓᵢ, ℓⱼ, γ, Aᵢ, Bᵢ, Pᵢ)
    sᵢ = 0

    for j in 0:floor(Int64, ((ℓᵢ + ℓⱼ) / 2))
        sᵢ +=
            cₖ((2 * j), ℓᵢ, ℓⱼ, (Pᵢ - Aᵢ), (Pᵢ - Bᵢ)) * doublefactorial(2 * j - 1) /
            (2 * γ)^j
    end
    sᵢ *= sqrt(π / γ)
    return sᵢ
end

function gaussianproduct(αᵢ, Rᵢ, αⱼ, Rⱼ, γ)
    P = [
        (αᵢ * Rᵢ[1] + αⱼ * Rⱼ[1]) / γ;
        (αᵢ * Rᵢ[2] + αⱼ * Rⱼ[2]) / γ;
        (αᵢ * Rᵢ[3] + αⱼ * Rⱼ[3]) / γ
    ]

    return P
end