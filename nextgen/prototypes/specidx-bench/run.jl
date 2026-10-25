# Spectral-index sampler benchmark: Commander3-style Gibbs/MH vs gradient-based samplers.
# Usage:  julia --project=. -t auto run.jl
# Env:    NPIX (default 96), NOISE (noise scale, default 1.0), SEED (default 1)

using Random, Statistics, Printf, MCMCDiagnosticTools
include(joinpath(@__DIR__, "src", "model.jl"))
include(joinpath(@__DIR__, "src", "samplers.jl"))

const NPIX  = parse(Int, get(ENV, "NPIX", "96"))
const NOISE = parse(Float64, get(ENV, "NOISE", "1.0"))
const SEED  = parse(Int, get(ENV, "SEED", "1"))
const PARAMS = ["a_cmb", "A_s", "A_d", "beta_s", "beta_d"]

const SAMPLERS = [
    ("MH chisq k=10",     (rng, p) -> gibbs_mh(rng, p; k = 10,  lnl = :chisq)),
    ("MH chisq k=100",    (rng, p) -> gibbs_mh(rng, p; k = 100, lnl = :chisq)),
    ("MH marginal k=10",  (rng, p) -> gibbs_mh(rng, p; k = 10,  lnl = :marginal)),
    ("Gibbs + NUTS",      (rng, p) -> gibbs_nuts(rng, p)),
    ("Joint NUTS",        (rng, p) -> joint_nuts(rng, p)),
    ("Collapsed NUTS",    (rng, p) -> collapsed_nuts(rng, p)),
]

essv(x) = ess(reshape(x, :, 1, 1))[1]

struct Res
    time::Float64
    ess::Vector{Float64}
    mean::Vector{Float64}
    sd::Vector{Float64}
    q::Matrix{Float64}        # 4 x 5 : 2.5, 16, 84, 97.5 percentiles
    n::Int
    nsolve::Int
    neval::Int
    ngrad::Int
end

function run_one(f, rng, p)
    t0 = time_ns()
    draws, cost = f(rng, p)
    dt = (time_ns() - t0) / 1e9
    e  = [essv(draws[:, j]) for j in 1:5]
    q  = reduce(hcat, [quantile(draws[:, j], [0.025, 0.16, 0.84, 0.975]) for j in 1:5])
    Res(dt, e, vec(mean(draws, dims = 1)), vec(std(draws, dims = 1)), q, size(draws, 1),
        cost.nsolve, cost.neval, cost.ngrad)
end

function main()
    pix = simulate(NPIX; seed = SEED, noise_scale = NOISE)
    @printf("Pixels: %d   noise scale: %.2f   threads: %d\n", NPIX, NOISE, Threads.nthreads())

    # compile everything once so timings exclude JIT
    for (_, f) in SAMPLERS
        run_one(f, Xoshiro(0), pix[end])
    end

    results = Dict{String,Vector{Res}}()
    for (si, (name, f)) in enumerate(SAMPLERS)
        rs = Vector{Res}(undef, NPIX)
        wall = @elapsed Threads.@threads :static for i in 1:NPIX
            rs[i] = run_one(f, Xoshiro(SEED * 1_000_003 + 10_007 * si + i), pix[i])
        end
        results[name] = rs
        @printf("  %-18s done in %6.1f s wall\n", name, wall)
    end

    ref = results["Collapsed NUTS"]
    tb  = [p.t for p in pix]
    bins = [(0.0, 1/3, "faint"), (1/3, 2/3, "mid"), (2/3, 1.0, "bright")]

    open(joinpath(@__DIR__, "summary_noise$(NOISE).txt"), "w") do io
        for out in (stdout, io)
            println(out, "\n=== Efficiency (medians over pixels; minESS = worst of the 5 parameters) ===")
            @printf(out, "%-18s %10s %10s %12s %14s %14s\n",
                    "sampler", "minESS/s", "ESS/1k it", "solves/ESS", "grad+eval/ESS", "time/pix [ms]")
            for (name, _) in SAMPLERS
                rs = results[name]
                m  = [minimum(r.ess) for r in rs]
                @printf(out, "%-18s %10.1f %10.1f %12.2f %14.1f %14.1f\n", name,
                        median(m ./ [r.time for r in rs]),
                        median(1000 .* m ./ [r.n for r in rs]),
                        median([r.nsolve for r in rs] ./ m),
                        median(([r.ngrad + r.neval for r in rs]) ./ m),
                        1000 * median([r.time for r in rs]))
            end

            println(out, "\n=== minESS per second by pixel brightness ===")
            @printf(out, "%-18s %10s %10s %10s\n", "sampler", bins[1][3], bins[2][3], bins[3][3])
            for (name, _) in SAMPLERS
                rs = results[name]
                vals = map(bins) do (lo, hi, _)
                    idx = findall(t -> lo <= t < hi || (hi == 1.0 && t == 1.0), tb)
                    median([minimum(rs[i].ess) / rs[i].time for i in idx])
                end
                @printf(out, "%-18s %10.1f %10.1f %10.1f\n", name, vals...)
            end

            println(out, "\n=== Which parameter limits each sampler (count of pixels) ===")
            for (name, _) in SAMPLERS
                worst = [PARAMS[argmin(r.ess)] for r in results[name]]
                counts = join(["$(k)=$(count(==(k), worst))" for k in PARAMS if count(==(k), worst) > 0], "  ")
                @printf(out, "%-18s %s\n", name, counts)
            end

            println(out, "\n=== Agreement with Collapsed NUTS (median |mean diff| / ref sd) ===")
            @printf(out, "%-18s %10s %10s %10s %10s\n", "sampler", "beta_s", "beta_d", "A_s", "A_d")
            for (name, _) in SAMPLERS
                rs = results[name]
                d(j) = median([abs(rs[i].mean[j] - ref[i].mean[j]) / ref[i].sd[j] for i in 1:NPIX])
                @printf(out, "%-18s %10.3f %10.3f %10.3f %10.3f\n", name, d(4), d(5), d(2), d(3))
            end

            println(out, "\n=== Coverage of the true values (nominal 68% / 95%) ===")
            for name in ("Collapsed NUTS", "Joint NUTS", "MH chisq k=10")
                rs = results[name]
                cov(j, lo, hi) = mean([rs[i].q[lo, j] <= pix[i].truth[j] <= rs[i].q[hi, j] for i in 1:NPIX])
                @printf(out, "%-18s beta_s %4.0f%% / %4.0f%%   beta_d %4.0f%% / %4.0f%%\n", name,
                        100cov(4, 2, 3), 100cov(4, 1, 4), 100cov(5, 2, 3), 100cov(5, 1, 4))
            end
        end
    end

    open(joinpath(@__DIR__, "results_noise$(NOISE).csv"), "w") do io
        println(io, "sampler,pixel,t,time_s,n,nsolve,neval,ngrad,",
                join(["ess_" * p for p in PARAMS], ","), ",", join(["mean_" * p for p in PARAMS], ","), ",",
                join(["sd_" * p for p in PARAMS], ","), ",", join(["true_" * p for p in PARAMS], ","))
        for (name, _) in SAMPLERS, (i, r) in enumerate(results[name])
            println(io, join(Any[name, i, pix[i].t, r.time, r.n, r.nsolve, r.neval, r.ngrad,
                                 r.ess..., r.mean..., r.sd..., pix[i].truth...], ","))
        end
    end
end

main()
