@inline function θ(l, lA, lB, PA, PB, r, g)
    θ = cₖ(l, lA, lB, PA, PB)
    θ *= factorial(l) * g^(r - l)
    θ /= factorial(r) * factorial(l - 2 * r)

    return θ
end

@inline function gi(l, lp, r, rp, i, lA, lB, Ai, Bi, Pi, gP, lC, lD, Ci, Di, Qi, gQ)
    δ = 1 / (4 * gP) + 1 / (4 * gQ)

    gi = (-1.0)^l
    gi *= θ(l, lA, lB, Pi - Ai, Pi - Bi, r, gP) * θ(lp, lC, lD, Qi - Ci, Qi - Di, rp, gQ)
    gi *= (-1.0)^i * (2 * δ)^(2 * (r + rp))
    gi *= factorial(l + lp - 2 * r - 2 * rp) * δ^i
    gi *= (Pi - Qi)^(l + lp - 2 * (r + rp + i))
    gi /= (4 * δ)^(l + lp) * factorial(i)
    gi /= factorial(l + lp - 2 * (r + rp + i))

    return gi
end

function Gxyz(lA, mA, nA, lB, mB, nB, lC, mC, nC, lD, mD, nD, a, b, c, d, RA, RB, RC, RD)
    gP = a + b
    gQ = c + d

    δ = 1 / (4 * gP) + 1 / (4 * gQ)

    RP = gaussianproduct(a, RA, b, RB, gP)
    RQ = gaussianproduct(c, RC, d, RD, gQ)

    AB = distance(RA, RB)
    CD = distance(RC, RD)
    PQ = distance(RP, RQ)

    Gxyz = 0.0

    for l = 0:(lA+lB)
        for r = 0:trunc(Int64, l / 2)
            for lp = 0:(lC+lD)
                for rp = 0:trunc(Int64, lp / 2)
                    for i = 0:trunc(Int64, (l + lp - 2 * r - 2 * rp) / 2)
                        gx = gi(l, lp, r, rp, i, lA, lB, RA[1], RB[1], RP[1], gP, lC, lD, RC[1], RD[1], RQ[1], gQ)
                        println("Calculating gx for l=$l, lp=$lp, r=$r, rp=$rp, i=$i...")
                        for m = 0:(mA+mB)
                            for s = 0:trunc(Int64, m / 2)
                                for mp = 0:(mC+mD)
                                    for sp = 0:trunc(Int64, mp / 2)
                                        for j = 0:trunc(Int64, (m + mp - 2 * s - 2 * sp) / 2)
                                            gy = gi(m, mp, s, sp, j, mA, mB, RA[2], RB[2], RP[2], gP, mC, mD, RC[2], RD[2], RQ[2], gQ)
                                            println("Calculating gy for m=$m, mp=$mp, s=$s, sp=$sp, j=$j...")
                                            for n = 0:(nA+nB)
                                                for t = 0:trunc(Int64, n / 2)
                                                    for np = 0:(nC+nD)
                                                        for tp = 0:trunc(Int64, np / 2)
                                                            for k = 0:trunc(Int64, (n + np - 2 * t - 2 * tp) /2)
                                                                gz = gi(n, np, t, tp, k, nA, nB, RA[3], RB[3], RP[3], gP, nC, nD, RC[3], RD[3], RQ[3], gQ)
                                                                println("Calculating gz for n=$n, np=$np, t=$t, tp=$tp, k=$k...")
                                                                ν = l + lp + m + mp + n + np - 2 * (r + rp + s + sp + t + tp) - (i + j + k)
                                                                F = boys(ν, PQ / (4 * δ))
                                                                Gxyz += gx * gy * gz * F
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    Gxyz *= (2 * π^2) / (gP * gQ)
    Gxyz *= sqrt(π / (gP + gQ))
    Gxyz *= exp(-(a * b * AB) / gP)
    Gxyz *= exp(-(c * d * CD) / gQ)

    return Gxyz
end

function repulsion(basis, molecule::Molecule)
    n = length(basis)
    G = zeros(n, n, n, n)

    for c in CartesianIndices(G)
        i, j, t, u = c[1], c[2], c[3], c[4]

        basisᵢ, basisⱼ, basisₜ, basisᵤ = basis[i], basis[j], basis[t], basis[u]
        Rᵢ, Rⱼ, Rₜ, Rᵤ = basisᵢ.R, basisⱼ.R, basisₜ.R, basisᵤ.R
        m, p, q, r = basisᵢ.size, basisⱼ.size, basisₜ.size, basisᵤ.size

        ℓᵢ, mᵢ, nᵢ = basisᵢ.ℓ, basisᵢ.m, basisᵢ.n
        ℓⱼ, mⱼ, nⱼ = basisⱼ.ℓ, basisⱼ.m, basisⱼ.n
        ℓₜ, mₜ, nₜ = basisₜ.ℓ, basisₜ.m, basisₜ.n
        ℓᵤ, mᵤ, nᵤ = basisᵤ.ℓ, basisᵤ.m, basisᵤ.n

        for e in CartesianIndices((m, p, q, r))
            k, l, s, v = e[1], e[2], e[3], e[4]

            αᵢ, αⱼ, αₜ, αᵤ = basisᵢ.α[k], basisⱼ.α[l], basisₜ.α[s], basisᵤ.α[v]
            dᵢ, dⱼ, dₜ, dᵤ = basisᵢ.d[k], basisⱼ.d[l], basisₜ.d[s], basisᵤ.d[v]
            Nᵢ, Nⱼ, Nₜ, Nᵤ = basisᵢ.N[k], basisⱼ.N[l], basisₜ.N[s], basisᵤ.N[v]

            println("Calculating G[$i, $j, $t, $u]...")
            tei = dᵢ * dⱼ * dₜ * dᵤ * Nᵢ * Nⱼ * Nₜ * Nᵤ
            tei *= Gxyz(ℓᵢ, mᵢ, nᵢ, ℓⱼ, mⱼ, nⱼ, ℓₜ, mₜ, nₜ, ℓᵤ, mᵤ, nᵤ, αᵢ, αⱼ, αₜ, αᵤ, Rᵢ, Rⱼ, Rₜ, Rᵤ)

            G[i, j, t, u] += tei
        end
    end

    return G
end