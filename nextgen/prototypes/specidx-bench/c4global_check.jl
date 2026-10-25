# Longer, multi-chain check of part B in c4test.jl
using Random, Statistics, Printf
redirect_stdout(devnull) do
    include(joinpath(@__DIR__, "c4test.jl"))   # defines functions; its main() output is discarded
end
for npix in (256, 1024)
    gp = simulate_global(npix, -3.05, 1.58; seed = 7 + npix)
    ex, _, _ = run_nuts(Xoshiro(1), GlobalTarget(gp, false), [PRIOR_S[1], PRIOR_D[1]]; nadapt = 500, nsamp = 6000)
    exm = reduce(vcat, permutedims.(ex))
    chains = Vector{Matrix{Float64}}(undef, 4)
    Threads.@threads :static for c in 1:4
        chains[c] = c4_global(Xoshiro(100 + c), gp; niter = 12000, nburn = 2000)
    end
    for (j, nm) in ((1, "beta_s"), (2, "beta_d"))
        m = [mean(ch[:, j]) for ch in chains]
        s = mean(std(ch[:, j]) for ch in chains)
        @printf("npix %5d %-6s exact %.4f ± %.4f | C4 chains mean %.4f, chain-to-chain se %.4f | bias %+.2f sd (±%.2f) | sd ratio %.2f\n",
                npix, nm, mean(exm[:, j]), std(exm[:, j]), mean(m), std(m) / 2,
                (mean(m) - mean(exm[:, j])) / std(exm[:, j]), std(m) / 2 / std(exm[:, j]), s / std(exm[:, j]))
    end
end
