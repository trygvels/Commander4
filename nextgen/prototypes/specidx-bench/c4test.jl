# Test of the Commander4 spectral-index move (compsep/mcmc.py + spectral_index.py):
# each Gibbs iteration draws amplitudes a ~ p(a | beta, d), then MCMC steps that
#   propose beta' (Gaussian random walk), redraw a' ~ p(a | beta', d) (constrained realization),
#   accept with  delta(-chi2/2) + delta(log prior),  revert beta and amplitudes on reject.
# Exact reference: amplitudes marginalized analytically ("collapsed"), sampled with NUTS.
# Also computed: the "ridge" (profile-likelihood) posterior, i.e. the collapsed posterior
# without the log det term, which the derivation in the README predicts this move targets.
#
# Usage: julia --project=. -t auto c4test.jl

using Random, Statistics, Printf, MCMCDiagnosticTools
include(joinpath(@__DIR__, "src", "model.jl"))
include(joinpath(@__DIR__, "src", "samplers.jl"))

essv(x) = ess(reshape(x, :, 1, 1))[1]

# ---------------- per-pixel spectral indices ----------------
function c4_mh(rng, p::Pixel; niter = 4000, nburn = 1000, numstep = 1)
    b = MVector(PRIOR_S[1], PRIOR_D[1])
    a = draw_amps(rng, p, b[1], b[2])
    step = MVector(0.05, 0.02)
    nacc, nprop, nsolve = 0, 0, 1
    out = Matrix{Float64}(undef, niter, 5)
    for it in 1:(nburn + niter)
        a = draw_amps(rng, p, b[1], b[2]); nsolve += 1          # amplitude group
        for _ in 1:numstep                                      # MCMC group
            cur = loglike(p, a, b[1], b[2]) + logprior(b[1], b[2])
            b2 = b .+ step .* randn(rng, 2)
            a2 = draw_amps(rng, p, b2[1], b2[2]); nsolve += 1   # re-solve with fluctuations
            new = loglike(p, a2, b2[1], b2[2]) + logprior(b2[1], b2[2])
            if log(rand(rng)) < new - cur
                b .= b2; a = a2; nacc += 1
            end
            nprop += 1
        end
        if it <= nburn && nprop >= 100                           # tuning (C4 itself uses a fixed sigma)
            acc = nacc / nprop
            step .*= acc > 0.5 ? 1.3 : (acc < 0.2 ? 0.6 : 1.0)
            nacc = 0; nprop = 0
        end
        it > nburn && (out[it - nburn, :] .= (a[1], a[2], a[3], b[1], b[2]))
    end
    out, nsolve
end

# ---------------- full-sky scalar spectral indices (what C4 samples today) ----------------
struct GlobalTarget; pix::Vector{Pixel}; ridge::Bool; end
function glogp(t::GlobalTarget, bs, bd)
    s = zero(promote_type(typeof(bs), typeof(bd)))
    for p in t.pix
        A  = mixing(bs, bd); WA = p.w .* A; M = A' * WA; bb = WA' * p.d
        dM = det(M)
        dM > 0 || return convert(typeof(dM), -Inf)
        s += 0.5 * dot(bb, M \ bb) - (t.ridge ? 0.0 : 0.5 * log(dM))
    end
    s + logprior(bs, bd)
end
LogDensityProblems.logdensity(t::GlobalTarget, x) = glogp(t, x[1], x[2])
LogDensityProblems.dimension(::GlobalTarget) = 2
LogDensityProblems.capabilities(::Type{GlobalTarget}) = LogDensityProblems.LogDensityOrder{0}()

function simulate_global(npix, bs, bd; seed = 7)
    rng = Xoshiro(seed)
    sig = SIGMA0; w = 1.0 ./ sig .^ 2
    map(1:npix) do i
        t = (i - 0.5) / npix
        a = SVector(70randn(rng), 10^(0.0 + 2.5t + 0.1randn(rng)), 10^(0.5 + 2.5t + 0.1randn(rng)))
        d = mixing(bs, bd) * a + sig .* randn(rng, SVector{NB,Float64})
        Pixel(d, w, SVector(a..., bs, bd), t)
    end
end

function c4_global(rng, pix; niter = 6000, nburn = 1500)
    b = MVector(PRIOR_S[1], PRIOR_D[1])
    amps(bv) = [draw_amps(rng, p, bv[1], bv[2]) for p in pix]
    ll(A, bv) = sum(loglike(pix[i], A[i], bv[1], bv[2]) for i in eachindex(pix)) + logprior(bv[1], bv[2])
    A = amps(b)
    step = MVector(0.02, 0.01); nacc = 0; nprop = 0
    out = Matrix{Float64}(undef, niter, 2)
    for it in 1:(nburn + niter)
        A = amps(b)
        cur = ll(A, b)
        b2 = b .+ step .* randn(rng, 2)
        A2 = amps(b2)
        new = ll(A2, b2)
        if log(rand(rng)) < new - cur
            b .= b2; A = A2; nacc += 1
        end
        nprop += 1
        if it <= nburn && nprop >= 50
            acc = nacc / nprop
            step .*= acc > 0.5 ? 1.3 : (acc < 0.2 ? 0.6 : 1.0)
            nacc = 0; nprop = 0
        end
        it > nburn && (out[it - nburn, :] .= b)
    end
    out
end

function main()
    println("=== A. Per-pixel indices: Commander4-style move vs exact (300 pixels, noise x1) ===")
    pix = simulate(300; seed = 1, noise_scale = 1.0)
    c4_mh(Xoshiro(0), pix[1]); collapsed_nuts(Xoshiro(0), pix[1])
    n = length(pix)
    c4 = Vector{Any}(undef, n); ref = Vector{Any}(undef, n)
    Threads.@threads :static for i in 1:n
        c4[i]  = c4_mh(Xoshiro(10 + i), pix[i])
        ref[i] = collapsed_nuts(Xoshiro(20 + i), pix[i])[1]
    end
    for (j, nm) in ((4, "beta_s"), (5, "beta_d"), (2, "A_s"), (3, "A_d"), (1, "a_cmb"))
        z  = [(mean(c4[i][1][:, j]) - mean(ref[i][:, j])) / std(ref[i][:, j]) for i in 1:n]
        sr = [std(c4[i][1][:, j]) / std(ref[i][:, j]) for i in 1:n]
        @printf("%-7s  mean shift (ref sd): median %+.3f  |median| %.3f   sd ratio C4/exact: median %.3f\n",
                nm, median(z), median(abs.(z)), median(sr))
    end
    ess4 = [minimum(essv(c4[i][1][:, j]) for j in 1:5) for i in 1:n]
    @printf("C4-style: minESS per 1k it %.1f, amplitude solves per ESS %.2f\n",
            median(1000 .* ess4 ./ 4000), median([c4[i][2] for i in 1:n] ./ ess4))
    bright = findall(p -> p.t > 2/3, pix)
    @printf("bright pixels only: beta_d mean shift median %+.3f ref sd\n",
            median([(mean(c4[i][1][:, 5]) - mean(ref[i][:, 5])) / std(ref[i][:, 5]) for i in bright]))

    println("\n=== B. Full-sky scalar indices (Commander4's current use): bias vs number of pixels ===")
    @printf("%6s %12s %12s %12s %12s %14s %14s\n", "npix", "exact b_s", "exact sd", "C4 b_s", "ridge b_s",
            "C4 bias [sd]", "C4 bias b_d [sd]")
    for npix in (16, 64, 256, 1024)
        gp = simulate_global(npix, -3.05, 1.58; seed = 7 + npix)
        ex, _, _ = run_nuts(Xoshiro(1), GlobalTarget(gp, false), [PRIOR_S[1], PRIOR_D[1]]; nadapt = 500, nsamp = 3000)
        rd, _, _ = run_nuts(Xoshiro(2), GlobalTarget(gp, true),  [PRIOR_S[1], PRIOR_D[1]]; nadapt = 500, nsamp = 3000)
        exm = reduce(vcat, permutedims.(ex)); rdm = reduce(vcat, permutedims.(rd))
        c = c4_global(Xoshiro(3), gp)
        e4 = essv(c[:, 1])
        @printf("%6d %12.4f %12.4f %12.4f %12.4f %14.2f %14.2f   (C4 ESS %.0f)\n", npix,
                mean(exm[:, 1]), std(exm[:, 1]), mean(c[:, 1]), mean(rdm[:, 1]),
                (mean(c[:, 1]) - mean(exm[:, 1])) / std(exm[:, 1]),
                (mean(c[:, 2]) - mean(exm[:, 2])) / std(exm[:, 2]), e4)
    end
end

main()
