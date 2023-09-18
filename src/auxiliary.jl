"""
This function calculate the double factorial of a number.
"""
function doublefactorial(number)
    fact = one(number)
    if number % 2 == 0
        for m = 1:number
            if m % 2 == 0
                fact *= m
            end
        end
    elseif number % 2 == 1
        for m = 1:number
            if m % 2 == 1
                fact *= m
            end
        end
    end

    return fact
end

"""Compute the square of the distance between two points.
            |A|² = (Aₐ - Bₐ)² - (Aₑ - Bₑ)² - (Aₒ - Bₒ)². 
The coefficients are multiplied by the (x,y,z) coordinates of the points.
"""
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