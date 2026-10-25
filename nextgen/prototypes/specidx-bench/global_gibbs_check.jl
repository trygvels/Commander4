# Independent check of the exact full-sky reference: Commander3-style Gibbs
# (amplitudes drawn given beta, then MH on beta given amplitudes). Shares no math with the
# collapsed formula, so agreement validates both.
using Random, Statistics, Printf
redirect_stdout(devnull) do
    include(joinpath(@__DIR__, "c4test.jl"))
end
function gibbs_global(rng, pix; niter = 40000, nburn = 4000, k = 20)
    b = MVector(PRIOR_S[1], PRIOR_D[1]); step = MVector(0.005, 0.002)
    out = Matrix{Float64}(undef, niter, 2); nacc = MVector(0, 0); nprop = MVector(0, 0)
    for it in 1:(nburn + niter)
        A = [draw_amps(rng, p, b[1], b[2]) for p in pix]
        ll(x, y) = sum(loglike(pix[i], A[i], x, y) for i in eachindex(pix)) + logprior(x, y)
        for j in 1:2
            cur = ll(b[1], b[2])
            for _ in 1:k
                b2 = copy(b); b2[j] += step[j] * randn(rng)
                new = ll(b2[1], b2[2])
                if log(rand(rng)) < new - cur
                    b .= b2; cur = new; nacc[j] += 1
                end
                nprop[j] += 1
            end
            if it <= nburn && nprop[j] >= 100
                acc = nacc[j] / nprop[j]; step[j] *= acc > 0.7 ? 1.5 : (acc < 0.3 ? 0.5 : 1.0)
                nacc[j] = 0; nprop[j] = 0
            end
        end
        it > nburn && (out[it - nburn, :] .= b)
    end
    out
end
for npix in (64, 256)
    gp = simulate_global(npix, -3.05, 1.58; seed = 7 + npix)
    ex, _, _ = run_nuts(Xoshiro(1), GlobalTarget(gp, false), [PRIOR_S[1], PRIOR_D[1]]; nadapt = 500, nsamp = 6000)
    exm = reduce(vcat, permutedims.(ex))
    ch = Vector{Matrix{Float64}}(undef, 2)
    Threads.@threads :static for c in 1:2
        ch[c] = gibbs_global(Xoshiro(500 + c), gp)
    end
    for (j, nm) in ((1, "beta_s"), (2, "beta_d"))
        allc = vcat(ch[1][:, j], ch[2][:, j]); e = essv(ch[1][:, j]) + essv(ch[2][:, j])
        @printf("npix %4d %-6s exact %.4f ± %.4f | Gibbs %.4f ± %.4f (ESS %.0f) | shift %+.2f sd (MC err ±%.2f) | sd ratio %.2f\n",
                npix, nm, mean(exm[:, j]), std(exm[:, j]), mean(allc), std(allc), e,
                (mean(allc) - mean(exm[:, j])) / std(exm[:, j]), std(allc) / sqrt(e) / std(exm[:, j]),
                std(allc) / std(exm[:, j]))
    end
end
